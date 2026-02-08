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

from ._c_backend import CBackendMixin
from .ContDynSys import ContDynSys
from ..utils.tools.chronometer import time_function_accumulate

if TYPE_CHECKING:
    from ..graphs.nx import SignedGraphNX as SignedGraph


class CoupledODEModel(CBackendMixin, ContDynSys):
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
        dynpath = getattr(sg, 'path_data', None)
        if dynpath is not None:
            dynpath = Path(dynpath) / "coupled_ode"
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
                    raise ValueError(f"Unknown coupling: {self.coupling_type!r}")
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
        self._check_c_backend_or_fallback()
        self.init_state(custom)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_state()
            if self.rndStr and not exName:
                exName = self.run_id
            self.sg._export_edgel_bin(exName=exName or self.run_id)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        if self.sini is None:
            self.init_ode_dynamics()

    # ------------------------------------------------------------------
    # C backend
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        coupling_str = self.coupling_type if isinstance(self.coupling_type, str) else "custom"
        local_str = self.local_type if isinstance(self.local_type, str) else "custom"
        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except (ValueError, AttributeError):
            datdir = getattr(self.sg, 'path_sgdata',
                             getattr(self.sg, 'path_data', Path.cwd()))
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
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
        return [getattr(self, 'sfout', None)]

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        verbose: bool = False,
        clean_export: bool = True,
        **kw: Any,
    ) -> None:
        self.check_attribute()
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.run_py(verbose=verbose)
