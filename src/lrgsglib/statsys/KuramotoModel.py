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

from ._c_backend import CBackendMixin
from .ContDynSys import ContDynSys
from ..utils.tools.chronometer import time_function_accumulate

if TYPE_CHECKING:
    from ..graphs.nx import SignedGraphNX as SignedGraph


class KuramotoModel(CBackendMixin, ContDynSys):
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

    # CBackendMixin configuration
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
        dynpath = getattr(sg, 'path_data', None)
        if dynpath is not None:
            dynpath = Path(dynpath) / "kuramoto"
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
        """Set natural frequencies."""
        if isinstance(omega, str):
            match omega:
                case "gaussian" | "normal":
                    self._omega = np.random.randn(self.N)
                case "uniform":
                    self._omega = np.random.uniform(-0.5, 0.5, self.N)
                case _:
                    raise ValueError(f"Unknown omega mode: {omega!r}")
        elif np.isscalar(omega):
            self._omega = np.full(self.N, float(omega), dtype=np.float64)
        else:
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
        diff = s[np.newaxis, :] - s[:, np.newaxis]  # diff[i, j] = theta_j - theta_i
        coupling_term = (self.coupling / self.N) * np.sum(A * np.sin(diff), axis=1)
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
    def init_kuramoto_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise phases and export data if C backend is selected."""
        self._check_c_backend_or_fallback()
        self.order_params = []
        self.init_state(custom)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_state()
            self.export_hfield()
            if self.rndStr and not exName:
                exName = self.run_id
            self.sg._export_edgel_bin(exName=exName or self.run_id)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        if self.sini is None:
            self.init_kuramoto_dynamics()

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
    # C backend integration
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        """Build argument list for KuramotoSimulator."""
        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except (ValueError, AttributeError):
            datdir = getattr(
                self.sg, 'path_sgdata',
                getattr(self.sg, 'path_data', Path.cwd())
            )
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
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
        # Read order parameter file if it exists
        op_path = self.dynpath / self.sg.get_p_fname('r', self.out_suffix)
        if op_path.exists():
            self.order_params = np.fromfile(op_path, dtype=np.float64).tolist()

    def _get_cleanup_paths(self) -> list[Path | None]:
        return [
            getattr(self, 'sfout', None),
            getattr(self, 'hfout', None),
        ]

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
        """Run Kuramoto dynamics.

        Parameters
        ----------
        verbose : bool
            Verbose output.
        clean_export : bool
            Remove exported files after C run.
        """
        self.check_attribute()
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.run_py(verbose=verbose)
