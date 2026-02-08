"""XY model (planar rotator) on signed graphs.

Each node *i* carries a phase theta_i in [0, 2*pi).  The Hamiltonian is::

    H = -sum_{(i,j)} A_ij * cos(theta_i - theta_j)

Updates use Metropolis-Hastings with a random angular perturbation.

Key observable: magnetisation m = |1/N * sum exp(i*theta)|, which
detects the Berezinskii-Kosterlitz-Thouless (BKT) transition in 2D.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from .._c_backend import CBackendMixin
from ..VecDynSys import VecDynSys
from ...utils.tools.chronometer import time_function_accumulate

if TYPE_CHECKING:
    from ...graphs.nx import SignedGraphNX as SignedGraph


class XYModel(CBackendMixin, VecDynSys):
    """XY (planar rotator) model with Metropolis dynamics.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    T : float
        Temperature.
    delta : float
        Maximum angular perturbation for proposals (in radians).
    steps : int
        Number of Monte Carlo sweeps.
    **kw
        Forwarded to :class:`VecDynSys`.

    Examples
    --------
    >>> xy = XYModel(lattice, T=0.5, steps=1000)
    >>> xy.init_xy_dynamics()
    >>> xy.run()
    >>> print(xy.magnetisation())
    """

    dyn_UVclass = "xy_model"

    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "XYSimulator{}"
    _allowed_c_keys: tuple[str, ...] = ("C0",)

    ene: list[float] = []
    magn: list[float] = []

    def __init__(
        self,
        sg: "SignedGraph",
        T: float = 1.0,
        delta: float = 0.5,
        *,
        steps: int | None = None,
        simref: float | None = None,
        save_observables: bool = False,
        **kw: Any,
    ) -> None:
        dynpath = getattr(sg, 'path_data', None)
        if dynpath is not None:
            dynpath = Path(dynpath) / "xy"
        super().__init__(
            sg,
            q=1,
            discrete=False,
            T=T,
            steps=steps,
            simref=simref,
            dynpath=dynpath,
            **kw,
        )
        self.delta = delta
        self.save_observables = save_observables
        self.ene = []
        self.magn = []
        self.sini: NDArray | None = None

    @property
    def _state_shape(self) -> tuple[int, ...]:
        return (self.N,)  # scalar angle per node

    # ------------------------------------------------------------------
    # Observables
    # ------------------------------------------------------------------
    def energy(self, s: NDArray | None = None) -> float:
        """Compute H = -sum A_ij cos(theta_i - theta_j)."""
        if s is None:
            s = self.s
        A = self._get_adj_matrix()
        E = 0.0
        for i in range(self.N):
            for j in range(i + 1, self.N):
                if A[i, j] != 0:
                    E -= A[i, j] * np.cos(s[i] - s[j])
        return E

    def magnetisation(self, s: NDArray | None = None) -> float:
        """Order parameter m = |1/N sum exp(i*theta)|."""
        if s is None:
            s = self.s
        return float(np.abs(np.mean(np.exp(1j * s))))

    # ------------------------------------------------------------------
    # Metropolis update
    # ------------------------------------------------------------------
    def ds1step(self, nd: int) -> None:
        """Single-node Metropolis update for XY model."""
        A = self._get_adj_matrix()
        current = self.s[nd]
        proposal = (current + np.random.uniform(-self.delta, self.delta)) % (2.0 * np.pi)

        # Energy change
        dE = 0.0
        for j in range(self.N):
            if A[nd, j] != 0:
                dE += A[nd, j] * (np.cos(current - self.s[j]) - np.cos(proposal - self.s[j]))

        if dE <= 0 or (self.T > 0 and np.random.uniform() < np.exp(-dE / self.T)):
            self.s[nd] = proposal

    # ------------------------------------------------------------------
    # Dynamics
    # ------------------------------------------------------------------
    def init_xy_dynamics(self, custom: Any = None, exName: str = "") -> None:
        self._check_c_backend_or_fallback()
        self.ene = []
        self.magn = []
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
            self.init_xy_dynamics()

    def run_py(self, verbose: bool = False, **kw: Any) -> None:
        """Monte Carlo sweeps with Metropolis."""
        for sweep in range(self.steps):
            nodes = np.random.permutation(self.N)
            for nd in nodes:
                self.ds1step(nd)
            if self.save_observables:
                self.ene.append(self.energy())
                self.magn.append(self.magnetisation())
            if self.savedyn:
                self.s_t.append(self.s.copy())

    # ------------------------------------------------------------------
    # C backend
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
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
            f"{self.T:.12g}",
            f"{self.delta:.12g}",
            str(datdir),
            syshape,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_suffix),
        ]

    def _get_cleanup_paths(self) -> list[Path | None]:
        return [getattr(self, 'sfout', None)]

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
