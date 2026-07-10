"""Reaction-diffusion dynamics on signed graphs.

Each node *i* carries a concentration u_i >= 0.  The dynamics are::

    du_i/dt = D * sum_j L_ij * u_j + f(u_i; params)

where *D* is the diffusion constant, *L* is the graph Laplacian, and
*f* is a local reaction function.

Built-in reaction types:

* ``"fisher_kpp"``: f(u) = r*u*(1-u)  (logistic growth)
* ``"bistable"``:   f(u) = u*(1-u)*(u - a) (Allee effect)
* ``"linear"``:     f(u) = r*u  (exponential growth/decay)
* ``"none"``:       f(u) = 0    (pure diffusion)
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from .._c_backend import CBackendMixin
from .._csr import build_graph_csr
from .._run_host import RunHostMixin
from ..ContDynSys import ContDynSys
from .defaults import (
    RD_DYN_SUBDIR,
    RD_INTEGRATOR_DEFAULT,
    RD_REACTION_DEFAULT,
    RD_SOLVER_NAME,
)

if TYPE_CHECKING:
    from ...graphs.nx import SignedGraphNX as SignedGraph


class ReactionDiffusionModel(RunHostMixin, CBackendMixin, ContDynSys):
    """Reaction-diffusion on a signed graph.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    D : float
        Diffusion constant.
    reaction : str
        Reaction type: ``"fisher_kpp"``, ``"bistable"``, ``"linear"``,
        ``"none"``.
    reaction_params : dict
        Parameters for the reaction function (e.g. ``{"r": 1.0}``).
    dt : float
        Integration time step.
    steps : int
        Number of integration steps.
    **kw
        Forwarded to :class:`ContDynSys`.

    Examples
    --------
    >>> rd = ReactionDiffusionModel(lattice, D=0.1, reaction="fisher_kpp",
    ...     reaction_params={"r": 1.0}, steps=1000)
    >>> rd.init_rd_dynamics()
    >>> rd.run()
    """

    dyn_UVclass = "reaction_diffusion"

    solver_name = RD_SOLVER_NAME

    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "ReactionDiffusionSimulator{}"
    _allowed_c_keys: tuple[str, ...] = ("C0",)

    def __init__(
        self,
        sg: "SignedGraph",
        D: float = 0.1,
        reaction: str = "fisher_kpp",
        reaction_params: dict[str, float] | None = None,
        *,
        dt: float = 0.01,
        integrator: str = "rk4",
        steps: int | None = None,
        simref: float | None = None,
        **kw: Any,
    ) -> None:
        # Graph-anchored dynamics tree (data/<graph>/reaction_diffusion/
        # N=...), like path_ising; flat <path_data>/<subdir> as a fallback.
        dynpath = getattr(sg, "path_reaction_diffusion", None)
        if dynpath is None:
            base = getattr(sg, "path_data", None)
            if base is not None:
                dynpath = Path(base) / RD_DYN_SUBDIR
        super().__init__(
            sg,
            dt=dt,
            integrator=integrator,
            steps=steps,
            simref=simref,
            dynpath=dynpath,
            **kw,
        )
        self.D = D
        self.reaction = reaction
        self.reaction_params = reaction_params or {}
        self._laplacian: NDArray | None = None
        self.sini: NDArray | None = None

    # ------------------------------------------------------------------
    # Laplacian
    # ------------------------------------------------------------------
    def _get_laplacian(self) -> NDArray:
        """Return the (cached) signed Laplacian matrix."""
        if self._laplacian is None:
            self._laplacian = self.sg.get_signed_laplacian().astype(np.float64)
        return self._laplacian

    # ------------------------------------------------------------------
    # Reaction functions
    # ------------------------------------------------------------------
    def _reaction(self, u: NDArray) -> NDArray:
        """Compute the local reaction term f(u)."""
        match self.reaction:
            case "fisher_kpp" | "fisher":
                r = self.reaction_params.get("r", 1.0)
                return r * u * (1.0 - u)
            case "bistable":
                a = self.reaction_params.get("a", 0.5)
                return u * (1.0 - u) * (u - a)
            case "linear":
                r = self.reaction_params.get("r", 1.0)
                return r * u
            case "none":
                return np.zeros_like(u)
            case _:
                raise ValueError(f"Unknown reaction type: {self.reaction!r}")

    # ------------------------------------------------------------------
    # ODE right-hand side
    # ------------------------------------------------------------------
    def rhs(self, s: NDArray, t: float) -> NDArray:
        """du/dt = -D * L @ u + f(u)."""
        L = self._get_laplacian()
        diffusion = -self.D * L @ s
        return diffusion + self._reaction(s)

    def _apply_constraints(self, s: NDArray) -> NDArray:
        """Clamp concentrations to [0, inf)."""
        return np.maximum(s, 0.0)

    # ------------------------------------------------------------------
    # State initialisation override
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise concentration field.

        Modes
        -----
        uniform / random:  U(0, 1).
        delta:             Single node at concentration 1.
        ones:              All nodes at 1.
        zero / zeros:      All nodes at 0.
        custom:            User-supplied array.
        """
        match self.ic:
            case "uniform" | "random" | "rand":
                self.s = np.random.uniform(0, 1, size=self.N)
            case "delta":
                self.s = np.zeros(self.N, dtype=np.float64)
                self.s[self.N // 2] = 1.0
            case "ones":
                self.s = np.ones(self.N, dtype=np.float64)
            case "zero" | "zeros":
                self.s = np.zeros(self.N, dtype=np.float64)
            case "custom":
                if custom is None:
                    raise ValueError("Provide a custom state array.")
                self.s = np.asarray(custom, dtype=np.float64).copy()
                if self.s.shape != (self.N,):
                    raise ValueError(
                        f"Shape mismatch: {self.s.shape} != ({self.N},)"
                    )
            case _:
                raise ValueError(f"Unknown IC: {self.ic!r}")

    # ------------------------------------------------------------------
    # Observables
    # ------------------------------------------------------------------
    def mean_concentration(self) -> float:
        return float(np.mean(self.s))

    def total_mass(self) -> float:
        return float(np.sum(self.s))

    # ------------------------------------------------------------------
    # Dynamics initialisation
    # ------------------------------------------------------------------
    def init_rd_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the concentration field.

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
            self.init_rd_dynamics()

    # ------------------------------------------------------------------
    # Run-dirname schema + per-run outputs (Phase C)
    # ------------------------------------------------------------------
    def _name_tokens(self) -> list:
        return [
            ("p", float(self.sg.pflip)),
            ("D", float(self.D)),
            *(
                (key, self.reaction_params[key])
                for key in sorted(self.reaction_params)
            ),
            ("dt", float(self.dt)),
            ("ns", int(self.steps)),
            ("rx", self.reaction, RD_REACTION_DEFAULT),
            ("intg", self.integrator, RD_INTEGRATOR_DEFAULT),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _cfg_model_block(self) -> dict[str, Any]:
        return {
            "class": type(self).__name__,
            "D": self.D,
            "reaction": self.reaction,
            "reaction_params": dict(self.reaction_params),
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
        # The C RD_DIR macro composes ``<datdir>/reaction_diffusion/
        # <syshape>/``; ``<syshapePth>/<run dirname>`` redirects the C
        # file exchange into the per-run directory.
        datdir = self._get_datdir_arg()
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        rundir = self._run_output_dir()
        if rundir is not None:
            syshape = f"{syshape}/{rundir.name}"
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.steps}",
            f"{self.D:.12g}",
            f"{self.dt:.12g}",
            self.reaction,
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
        """Integrate via the in-process ``_rd_native`` RK4 kernel.

        Marshals the graph into the engine-agnostic CSR triple (NX + GT) and
        runs the same ``rd_rhs`` + ``rk4_step`` kernels the C subprocess uses,
        with the identical per-step clamp to [0, inf). Deterministic, so it
        matches the Python integrator numerically. Mirrors
        ``KuramotoModel._run_pybind``.
        """
        from .ccore import _rd_native

        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        u0 = np.ascontiguousarray(self.s, dtype=np.float64)
        # Resolve reaction params exactly as the Python ``_reaction`` does.
        r = float(self.reaction_params.get("r", 1.0))
        a = float(self.reaction_params.get("a", 0.5))
        save_mean = bool(getattr(self, "savedyn", False))
        u, mean = _rd_native.rd_sampling(
            u0,
            ni,
            nw,
            nptr,
            float(self.D),
            float(self.dt),
            str(self.reaction),
            r,
            a,
            int(self.steps),
            save_mean,
        )
        self.s = u
        if save_mean:
            self.mean_concentrations = mean.tolist()

    # Public interface: ``run()`` is inherited from RunHostMixin (backend
    # resolution, solver dispatch, output lifecycle, C-export cleanup).


# Register ReactionDiffusionModel's solver backends (py/pb/C) in the registry.
from . import _solvers  # noqa: E402,F401
