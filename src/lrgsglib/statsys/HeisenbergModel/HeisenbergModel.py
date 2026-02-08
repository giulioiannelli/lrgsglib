"""3D Heisenberg model on signed graphs.

Each node *i* carries a unit vector **n**_i = (nx, ny, nz) on the
2-sphere.  The Hamiltonian is::

    H = -sum_{(i,j)} A_ij * (n_i . n_j)

Updates use Metropolis-Hastings with random rotations.

By the Mermin-Wagner theorem, there is no spontaneous magnetisation
in 2D at any T > 0.
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


class HeisenbergModel(CBackendMixin, VecDynSys):
    """3D Heisenberg model with Metropolis dynamics.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    T : float
        Temperature.
    delta : float
        Maximum rotation angle for proposals.
    steps : int
        Number of Monte Carlo sweeps.
    **kw
        Forwarded to :class:`VecDynSys`.

    Examples
    --------
    >>> hberg = HeisenbergModel(lattice, T=0.5, steps=1000)
    >>> hberg.init_heisenberg_dynamics()
    >>> hberg.run()
    >>> print(hberg.magnetisation())
    """

    dyn_UVclass = "heisenberg_model"

    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "HeisenbergSimulator{}"
    _allowed_c_keys: tuple[str, ...] = ("C0",)

    ene: list[float] = []
    magn: list[float] = []

    def __init__(
        self,
        sg: "SignedGraph",
        T: float = 1.0,
        delta: float = 0.3,
        *,
        steps: int | None = None,
        simref: float | None = None,
        save_observables: bool = False,
        **kw: Any,
    ) -> None:
        dynpath = getattr(sg, 'path_data', None)
        if dynpath is not None:
            dynpath = Path(dynpath) / "heisenberg"
        super().__init__(
            sg,
            q=3,
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
        return (self.N, 3)  # 3D unit vector per node

    # ------------------------------------------------------------------
    # State initialisation override
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise random unit vectors on S^2."""
        match self.ic:
            case "uniform" | "random" | "rand":
                # Random points on the sphere (Marsaglia)
                vecs = np.random.randn(self.N, 3)
                norms = np.linalg.norm(vecs, axis=1, keepdims=True)
                norms = np.where(norms > 0, norms, 1.0)
                self.s = (vecs / norms).astype(np.float64)
            case "zero" | "zeros":
                # All pointing along z
                self.s = np.zeros((self.N, 3), dtype=np.float64)
                self.s[:, 2] = 1.0
            case "custom":
                if custom is None:
                    raise ValueError("Provide a custom state array.")
                self.s = np.asarray(custom, dtype=np.float64).copy()
                if self.s.shape != (self.N, 3):
                    raise ValueError(f"Shape mismatch: {self.s.shape} != ({self.N}, 3)")
                # Normalise
                norms = np.linalg.norm(self.s, axis=1, keepdims=True)
                self.s /= np.where(norms > 0, norms, 1.0)
            case _:
                raise ValueError(f"Unknown IC: {self.ic!r}")

    # ------------------------------------------------------------------
    # Observables
    # ------------------------------------------------------------------
    def energy(self, s: NDArray | None = None) -> float:
        """Compute H = -sum A_ij (n_i . n_j)."""
        if s is None:
            s = self.s
        A = self._get_adj_matrix()
        E = 0.0
        for i in range(self.N):
            for j in range(i + 1, self.N):
                if A[i, j] != 0:
                    E -= A[i, j] * np.dot(s[i], s[j])
        return E

    def magnetisation(self, s: NDArray | None = None) -> float:
        """Magnitude of mean spin: |1/N sum n_i|."""
        if s is None:
            s = self.s
        m_vec = np.mean(s, axis=0)
        return float(np.linalg.norm(m_vec))

    def magnetisation_vector(self, s: NDArray | None = None) -> NDArray:
        """Mean spin vector."""
        if s is None:
            s = self.s
        return np.mean(s, axis=0)

    # ------------------------------------------------------------------
    # Metropolis update
    # ------------------------------------------------------------------
    def _propose_rotation(self, v: NDArray) -> NDArray:
        """Propose a small random rotation of unit vector v."""
        perturb = v + self.delta * np.random.uniform(-1, 1, size=3)
        norm = np.linalg.norm(perturb)
        if norm > 0:
            return perturb / norm
        return v.copy()

    def ds1step(self, nd: int) -> None:
        """Single-node Metropolis update for Heisenberg."""
        A = self._get_adj_matrix()
        current = self.s[nd].copy()
        proposal = self._propose_rotation(current)

        # Energy change
        dE = 0.0
        for j in range(self.N):
            if A[nd, j] != 0:
                dE += A[nd, j] * (np.dot(current, self.s[j]) - np.dot(proposal, self.s[j]))

        if dE <= 0 or (self.T > 0 and np.random.uniform() < np.exp(-dE / self.T)):
            self.s[nd] = proposal

    # ------------------------------------------------------------------
    # Dynamics
    # ------------------------------------------------------------------
    def init_heisenberg_dynamics(self, custom: Any = None, exName: str = "") -> None:
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
            self.init_heisenberg_dynamics()

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
