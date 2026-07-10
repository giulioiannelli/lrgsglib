"""Generic coupled ODE dynamics on signed graphs.

Each node *i* carries a scalar state x_i.  The dynamics are::

    dx_i/dt = f_local(x_i) + sum_j A_ij * g_coupling(x_i, x_j)

where ``f_local`` is a per-node autonomous term and ``g_coupling`` is
a pairwise coupling function mediated by the adjacency matrix.

Built-in coupling functions:

* ``"linear"``:     g(x_i, x_j) = x_j - x_i  (diffusive)
* ``"product"``:    g(x_i, x_j) = x_j * x_i

Built-in local functions:

* ``"none"``:       f(x) = 0
* ``"linear"``:     f(x) = a*x
* ``"logistic"``:   f(x) = r*x*(1 - x/K)
* ``"fitzhugh"``:   f(x) = x - x^3/3  (FitzHugh-Nagumo)

Custom functions can be provided as callables (Python backend only).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
from numpy.typing import NDArray

from .._c_backend import CBackendMixin
from .._csr import build_graph_csr
from .._run_host import RunHostMixin
from ..ContDynSys import ContDynSys
from .defaults import (
    CODE_COUPLING_DEFAULT,
    CODE_DYN_SUBDIR,
    CODE_INTEGRATOR_DEFAULT,
    CODE_LOCAL_DEFAULT,
    CODE_SOLVER_NAME,
)

if TYPE_CHECKING:
    from ...graphs.nx import SignedGraphNX as SignedGraph


def _render_callable(value: Any) -> Any:
    """Resolve a coupling/local spec for naming: callables become
    ``"custom"``-style labels, strings pass through."""
    if callable(value):
        return f"custom:{getattr(value, '__name__', 'callable')}"
    return value


class CoupledODEModel(RunHostMixin, CBackendMixin, ContDynSys):
    """Generic coupled ODE system on a signed graph.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    coupling_type : str or callable
        Coupling function type or callable ``g(x_i, x_j) -> float``.
    local_type : str or callable
        Local function type or callable ``f(x_i) -> float``.
    coupling_strength : float
        Multiplier for the coupling term.
    local_params : dict
        Parameters for built-in local functions.
    dt : float
        Integration time step.
    steps : int
        Number of integration steps.
    **kw
        Forwarded to :class:`ContDynSys`.

    Examples
    --------
    >>> ode = CoupledODEModel(lattice, coupling_type="linear",
    ...     local_type="fitzhugh", steps=500)
    >>> ode.init_ode_dynamics()
    >>> ode.run()
    """

    dyn_UVclass = "coupled_ode"

    solver_name = CODE_SOLVER_NAME

    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "CoupledODESimulator{}"
    _allowed_c_keys: tuple[str, ...] = ("C0",)

    def __init__(
        self,
        sg: "SignedGraph",
        coupling_type: str | Callable = "linear",
        local_type: str | Callable = "none",
        coupling_strength: float = 1.0,
        local_params: dict[str, float] | None = None,
        *,
        dt: float = 0.01,
        integrator: str = "rk4",
        steps: int | None = None,
        simref: float | None = None,
        **kw: Any,
    ) -> None:
        # Graph-anchored dynamics tree (data/<graph>/coupled_ode/N=...),
        # like path_ising; flat <path_data>/coupled_ode as a fallback.
        dynpath = getattr(sg, "path_coupled_ode", None)
        if dynpath is None:
            base = getattr(sg, "path_data", None)
            if base is not None:
                dynpath = Path(base) / CODE_DYN_SUBDIR
        super().__init__(
            sg,
            dt=dt,
            integrator=integrator,
            steps=steps,
            simref=simref,
            dynpath=dynpath,
            **kw,
        )
        self.coupling_type = coupling_type
        self.local_type = local_type
        self.coupling_strength = coupling_strength
        self.local_params = local_params or {}
        self.sini: NDArray | None = None

    # ------------------------------------------------------------------
    # Coupling functions
    # ------------------------------------------------------------------
    def _g_coupling(self, xi: float, xj: float) -> float:
        """Coupling function g(x_i, x_j)."""
        if callable(self.coupling_type):
            return self.coupling_type(xi, xj)
        match self.coupling_type:
            case "linear" | "diffusive":
                return xj - xi
            case "product":
                return xj * xi
            case _:
                raise ValueError(f"Unknown coupling: {self.coupling_type!r}")

    def _g_coupling_vec(self, s: NDArray) -> NDArray:
        """Vectorised coupling: for each node i, sum_j A_ij * g(x_i, x_j)."""
        A = self._get_adj_matrix()
        N = self.N
        result = np.zeros(N, dtype=np.float64)
        if callable(self.coupling_type):
            for i in range(N):
                for j in range(N):
                    if A[i, j] != 0:
                        result[i] += A[i, j] * self.coupling_type(s[i], s[j])
        else:
            match self.coupling_type:
                case "linear" | "diffusive":
                    # sum_j A_ij (x_j - x_i) = A @ x - deg * x_i
                    result = A @ s - np.abs(A).sum(axis=1) * s
                case "product":
                    result = (A * s[np.newaxis, :]) @ s
                case _:
                    raise ValueError(
                        f"Unknown coupling: {self.coupling_type!r}"
                    )
        return self.coupling_strength * result

    # ------------------------------------------------------------------
    # Local functions
    # ------------------------------------------------------------------
    def _f_local(self, s: NDArray) -> NDArray:
        """Local function f(x_i) applied elementwise."""
        if callable(self.local_type):
            return np.array([self.local_type(xi) for xi in s])
        match self.local_type:
            case "none":
                return np.zeros_like(s)
            case "linear":
                a = self.local_params.get("a", 1.0)
                return a * s
            case "logistic":
                r = self.local_params.get("r", 1.0)
                K = self.local_params.get("K", 1.0)
                return r * s * (1.0 - s / K)
            case "fitzhugh" | "fhn":
                return s - s**3 / 3.0
            case _:
                raise ValueError(f"Unknown local type: {self.local_type!r}")

    # ------------------------------------------------------------------
    # ODE right-hand side
    # ------------------------------------------------------------------
    def rhs(self, s: NDArray, t: float) -> NDArray:
        return self._f_local(s) + self._g_coupling_vec(s)

    # ------------------------------------------------------------------
    # Dynamics initialisation
    # ------------------------------------------------------------------
    def init_ode_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the node states.

        The C-backend file exchange (state / edgelist exports) now
        happens inside ``run()`` (``_begin_outputs``), into the per-run
        output directory; ``exName`` is kept for backward compatibility
        and ignored (exports are named by ``run_id``).
        """
        self._check_c_backend_or_fallback()
        self.init_state(custom)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        if self.sini is None:
            self.init_ode_dynamics()

    # ------------------------------------------------------------------
    # Run-dirname schema + per-run outputs (Phase C)
    # ------------------------------------------------------------------
    def _name_tokens(self) -> list:
        cpl = "custom" if callable(self.coupling_type) else self.coupling_type
        loc = "custom" if callable(self.local_type) else self.local_type
        return [
            ("p", float(self.sg.pflip)),
            ("K", float(self.coupling_strength)),
            *(
                (key, self.local_params[key])
                for key in sorted(self.local_params)
            ),
            ("dt", float(self.dt)),
            ("ns", int(self.steps)),
            ("cpl", cpl, CODE_COUPLING_DEFAULT),
            ("loc", loc, CODE_LOCAL_DEFAULT),
            ("intg", self.integrator, CODE_INTEGRATOR_DEFAULT),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _cfg_model_block(self) -> dict[str, Any]:
        return {
            "class": type(self).__name__,
            "coupling_type": _render_callable(self.coupling_type),
            "local_type": _render_callable(self.local_type),
            "coupling_strength": self.coupling_strength,
            "local_params": dict(self.local_params),
            "dt": self.dt,
            "integrator": self.integrator,
            "steps": self.steps,
            "ic": self.ic,
        }

    def _begin_outputs(self) -> None:
        self._open_run_outputs()

    def _persist_observables(self) -> None:
        self._flush_run_outputs()

    # ------------------------------------------------------------------
    # C backend
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        # The C CODE_DIR macro composes ``<datdir>/coupled_ode/<syshape>/``;
        # ``<syshapePth>/<run dirname>`` redirects the C file exchange
        # into the per-run directory.
        coupling_str = (
            self.coupling_type
            if isinstance(self.coupling_type, str)
            else "custom"
        )
        local_str = (
            self.local_type if isinstance(self.local_type, str) else "custom"
        )
        datdir = self._get_datdir_arg()
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        rundir = self._run_output_dir()
        if rundir is not None:
            syshape = f"{syshape}/{rundir.name}"
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.steps}",
            f"{self.coupling_strength:.12g}",
            f"{self.dt:.12g}",
            coupling_str,
            local_str,
            str(datdir),
            syshape,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_suffix),
        ]

    def _get_cleanup_paths(self) -> list[Path | None]:
        return [getattr(self, "sfout", None)]

    # ------------------------------------------------------------------
    # Native (pybind11) backend
    # ------------------------------------------------------------------
    def _run_pybind(self) -> None:
        """Integrate via the in-process ``_code_native`` RK4 kernel.

        Marshals the graph into the engine-agnostic CSR triple (NX + GT) and runs
        the same ``code_rhs`` / ``rk4_step`` kernel the C subprocess uses.
        Deterministic, so it matches the Python integrator numerically. Mirrors
        ``KuramotoModel._run_pybind``. Callable coupling/local functions are
        not supported by the native kernel (string types only).
        """
        from .ccore import _code_native

        if callable(self.coupling_type) or callable(self.local_type):
            raise TypeError(
                "_code_native supports only string coupling/local types; "
                "use runlang='py' for callable functions."
            )

        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        x0 = np.ascontiguousarray(self.s, dtype=np.float64)
        # Local-function parameters: defaults mirror the C kernel hardcoded
        # values, overridable via ``local_params`` (matching ``_f_local``).
        a = float(self.local_params.get("a", 1.0))
        r = float(self.local_params.get("r", 1.0))
        K = float(self.local_params.get("K", 1.0))
        x = _code_native.code_sampling(
            x0,
            ni,
            nw,
            nptr,
            float(self.coupling_strength),
            str(self.coupling_type),
            str(self.local_type),
            a,
            r,
            K,
            float(self.dt),
            int(self.steps),
        )
        self.s = x

    # Public interface: ``run()`` is inherited from RunHostMixin (backend
    # resolution, solver dispatch, output lifecycle, C-export cleanup).


# Register CoupledODEModel's solver backends (py/pb/C) in the shared registry.
from . import _solvers  # noqa: E402,F401
