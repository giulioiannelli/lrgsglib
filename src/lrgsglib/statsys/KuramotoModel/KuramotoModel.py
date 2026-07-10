"""Kuramoto model of coupled phase oscillators on signed graphs.

Each node *i* carries a phase theta_i in [0, 2*pi). The dynamics are::

    d(theta_i)/dt = omega_i + (K/N) * sum_j A_ij * sin(theta_j - theta_i)

where *K* is the global coupling strength, *omega_i* are natural
frequencies, and *A_ij* is the signed adjacency matrix.

Key observable: the complex order parameter
``r = |1/N * sum_j exp(i * theta_j)|``, which equals 1 for full
synchronisation and ~0 for incoherence.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from ...config.const import BIN
from .._c_backend import CBackendMixin
from .._csr import build_graph_csr
from .._run_host import RunHostMixin
from ..ContDynSys import ContDynSys
from .defaults import (
    KURAMOTO_DYN_SUBDIR,
    KURAMOTO_INTEGRATOR_DEFAULT,
    KURAMOTO_ORDER_FBASE,
    KURAMOTO_SOLVER_NAME,
)

if TYPE_CHECKING:
    from ...graphs.nx import SignedGraphNX as SignedGraph


class KuramotoModel(RunHostMixin, CBackendMixin, ContDynSys):
    """Kuramoto coupled oscillators on a signed graph.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    coupling : float
        Global coupling strength *K*.
    omega : float, NDArray, or str
        Natural frequencies.  A scalar sets all nodes equal; an array
        sets per-node frequencies; ``"gaussian"`` samples from N(0, 1);
        ``"uniform"`` samples from U(-0.5, 0.5).
    dt : float
        Integration time step.
    integrator : {"euler", "rk4", "rk45"}
        ODE integration method.
    steps : int
        Number of integration steps.
    **kw
        Forwarded to :class:`ContDynSys`.

    Examples
    --------
    >>> km = KuramotoModel(lattice, coupling=2.0, steps=1000)
    >>> km.init_kuramoto_dynamics()
    >>> km.run()
    >>> print(km.order_parameter())
    """

    dyn_UVclass = "kuramoto"

    solver_name = KURAMOTO_SOLVER_NAME

    # CBackendMixin configuration
    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "KuramotoSimulator{}"
    _allowed_c_keys: tuple[str, ...] = ("C0", "C1")

    # Observable storage
    order_params: list[float] = []

    def __init__(
        self,
        sg: "SignedGraph",
        coupling: float = 1.0,
        omega: float | NDArray | str = 0.0,
        *,
        dt: float = 0.01,
        integrator: str = "rk4",
        steps: int | None = None,
        simref: float | None = None,
        save_order_param: bool = False,
        **kw: Any,
    ) -> None:
        # Graph-anchored dynamics tree (data/<graph>/kuramoto/N=...),
        # like path_ising; flat <path_data>/kuramoto as a fallback.
        dynpath = getattr(sg, "path_kuramoto", None)
        if dynpath is None:
            base = getattr(sg, "path_data", None)
            if base is not None:
                dynpath = Path(base) / KURAMOTO_DYN_SUBDIR
        super().__init__(
            sg,
            dt=dt,
            integrator=integrator,
            steps=steps,
            simref=simref,
            dynpath=dynpath,
            **kw,
        )
        self.coupling = coupling
        self._init_omega(omega)
        self.save_order_param = save_order_param
        self.order_params = []
        self.sini: NDArray | None = None

    def _init_omega(self, omega: float | NDArray | str) -> None:
        """Set natural frequencies; remember the constructor's spec for
        the ``om=`` run-dirname token (scalar float / mode string /
        ``"custom"`` for an explicit array)."""
        if isinstance(omega, str):
            self._omega_spec: float | str = omega
            match omega:
                case "gaussian" | "normal":
                    self._omega = np.random.randn(self.N)
                case "uniform":
                    self._omega = np.random.uniform(-0.5, 0.5, self.N)
                case _:
                    raise ValueError(f"Unknown omega mode: {omega!r}")
        elif np.isscalar(omega):
            self._omega_spec = float(omega)
            self._omega = np.full(self.N, float(omega), dtype=np.float64)
        else:
            self._omega_spec = "custom"
            arr = np.asarray(omega, dtype=np.float64)
            if arr.shape != (self.N,):
                raise ValueError(
                    f"omega array shape {arr.shape} doesn't match N={self.N}"
                )
            self._omega = arr.copy()

    @property
    def omega(self) -> NDArray:
        """Natural frequencies (read-only)."""
        return self._omega

    # ------------------------------------------------------------------
    # ODE right-hand side
    # ------------------------------------------------------------------
    def rhs(self, s: NDArray, t: float) -> NDArray:
        """Kuramoto ODE: dtheta_i/dt = omega_i + (K/N) sum_j A_ij sin(theta_j - theta_i)."""
        A = self._get_adj_matrix()
        # sin(theta_j - theta_i) for all pairs, vectorised
        diff = (
            s[np.newaxis, :] - s[:, np.newaxis]
        )  # diff[i, j] = theta_j - theta_i
        coupling_term = (self.coupling / self.N) * np.sum(
            A * np.sin(diff), axis=1
        )
        return self._omega + coupling_term

    def _apply_constraints(self, s: NDArray) -> NDArray:
        """Wrap phases to [0, 2*pi)."""
        return s % (2.0 * np.pi)

    # ------------------------------------------------------------------
    # Observables
    # ------------------------------------------------------------------
    def order_parameter(self, s: NDArray | None = None) -> float:
        """Compute the Kuramoto order parameter r = |1/N sum exp(i*theta)|."""
        if s is None:
            s = self.s
        return float(np.abs(np.mean(np.exp(1j * s))))

    def mean_phase(self, s: NDArray | None = None) -> float:
        """Mean phase psi = arg(1/N sum exp(i*theta))."""
        if s is None:
            s = self.s
        return float(np.angle(np.mean(np.exp(1j * s))))

    # ------------------------------------------------------------------
    # Dynamics initialisation
    # ------------------------------------------------------------------
    def init_kuramoto_dynamics(
        self, custom: Any = None, exName: str = ""
    ) -> None:
        """Initialise phases.

        The C-backend file exchange (state / hfield / edgelist exports)
        now happens inside ``run()`` (``_begin_outputs``), into the
        per-run output directory; ``exName`` is kept for backward
        compatibility and ignored (exports are named by ``run_id``).
        """
        self._check_c_backend_or_fallback()
        self.order_params = []
        self.init_state(custom)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        if self.sini is None:
            self.init_kuramoto_dynamics()

    # ------------------------------------------------------------------
    # Run-dirname schema + per-run outputs (Phase C)
    # ------------------------------------------------------------------
    def _name_tokens(self) -> list:
        return [
            ("p", float(self.sg.pflip)),
            ("K", float(self.coupling)),
            ("om", self._omega_spec),
            ("dt", float(self.dt)),
            ("ns", int(self.steps)),
            ("intg", self.integrator, KURAMOTO_INTEGRATOR_DEFAULT),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _cfg_model_block(self) -> dict[str, Any]:
        return {
            "class": type(self).__name__,
            "coupling": self.coupling,
            "omega": self._omega_spec,
            "dt": self.dt,
            "integrator": self.integrator,
            "steps": self.steps,
            "save_order_param": self.save_order_param,
            "ic": self.ic,
        }

    def _begin_outputs(self) -> None:
        self.order_params = []
        self._open_run_outputs()

    def _export_c_inputs(self) -> None:
        super()._export_c_inputs()
        # KuramotoSimulator reads the natural-frequency file from the
        # hfield slot (the legacy convention: ``self.field`` is exported).
        self.export_hfield()

    def _persist_observables(self) -> None:
        if self.savedisk and self.save_order_param and len(self.order_params):
            rundir = self._run_output_dir()
            rundir.mkdir(parents=True, exist_ok=True)
            np.asarray(self.order_params, dtype=np.float64).tofile(
                rundir / f"{KURAMOTO_ORDER_FBASE}{BIN}"
            )
        self._flush_run_outputs()

    # ------------------------------------------------------------------
    # Python backend override (adds order parameter recording)
    # ------------------------------------------------------------------
    def run_py(self, verbose: bool = False, **kw: Any) -> None:
        """Integrate with order parameter tracking."""
        s = self.s.copy()
        for step in range(self.steps):
            s = self._integrate_step(s, step * self.dt)
            s = self._apply_constraints(s)
            if self.savedyn:
                self.s_t.append(s.copy())
            if self.save_order_param:
                self.order_params.append(self.order_parameter(s))
        self.s = s

    # ------------------------------------------------------------------
    # Native (pybind11) backend
    # ------------------------------------------------------------------
    def _run_pybind(self) -> None:
        """Integrate via the in-process ``_kuramoto_native`` RK4 kernel.

        Marshals the graph into the engine-agnostic CSR triple (NX + GT) and runs
        the same ``kuramoto_step_rk4`` kernel the C subprocess uses. Deterministic,
        so it matches the Python integrator numerically. Mirrors
        ``VoterModel._run_pybind``.
        """
        from .ccore import _kuramoto_native

        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        theta0 = np.ascontiguousarray(self.s, dtype=np.float64)
        omega = np.ascontiguousarray(self._omega, dtype=np.float64)
        theta, order = _kuramoto_native.kuramoto_sampling(
            theta0,
            ni,
            nw,
            nptr,
            omega,
            float(self.coupling),
            float(self.dt),
            int(self.steps),
            bool(self.save_order_param),
        )
        self.s = theta
        if self.save_order_param:
            self.order_params = order.tolist()

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        """Build argument list for KuramotoSimulator.

        The C ``KUR_DIR`` macro composes ``<datdir>/kuramoto/<syshape>/``;
        passing ``<syshapePth>/<run dirname>`` as the syshape argument
        redirects the whole C file exchange into the per-run directory.
        """
        datdir = self._get_datdir_arg()
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        rundir = self._run_output_dir()
        if rundir is not None:
            syshape = f"{syshape}/{rundir.name}"
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.steps}",
            f"{self.coupling:.12g}",
            f"{self.dt:.12g}",
            str(datdir),
            syshape,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_suffix),
        ]

    def run_cprogram(self, verbose: bool = False) -> None:
        """Execute C backend and read order parameter output."""
        super().run_cprogram(verbose)
        # Read the C-written order-parameter file from the per-run dir.
        op_path = self._c_exchange_dir() / self.sg.get_p_fname(
            "r", self.out_suffix
        )
        if op_path.exists():
            self.order_params = np.fromfile(op_path, dtype=np.float64).tolist()

    def _get_cleanup_paths(self) -> list[Path | None]:
        return [
            getattr(self, "sfout", None),
            getattr(self, "hfout", None),
        ]

    # Public interface: ``run()`` is inherited from RunHostMixin (backend
    # resolution, solver dispatch, output lifecycle, C-export cleanup).


# Register KuramotoModel's solver backends (py/pb/C) in the shared registry.
from . import _solvers  # noqa: E402,F401
