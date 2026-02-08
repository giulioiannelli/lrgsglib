"""Multi-species competition model on signed graphs.

Each node *i* carries a vector state sigma_i = (sigma_i^1, ..., sigma_i^k)
where each component is in {0, 1, ..., q-1}.  This generalises
multi-species contact/competition processes.

The Hamiltonian extends the Potts model to k coupled species::

    H = -sum_{(i,j)} sum_s interaction_matrix[s1, s2] * delta(sigma_i^s1, sigma_j^s2)

Updates use Metropolis-Hastings, flipping one component per proposal.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from ._c_backend import CBackendMixin
from .VecDynSys import VecDynSys
from ..utils.tools.chronometer import time_function_accumulate

if TYPE_CHECKING:
    from ..graphs.nx import SignedGraphNX as SignedGraph


class MultiSpeciesModel(CBackendMixin, VecDynSys):
    """Multi-species competition model with Metropolis dynamics.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    species : int
        Number of species (k).
    q_per_species : int or list[int]
        Number of states per species component.
    T : float
        Temperature.
    interaction_matrix : NDArray, optional
        Species interaction matrix (k x k).  Defaults to identity.
    steps : int
        Number of Monte Carlo sweeps.
    **kw
        Forwarded to :class:`VecDynSys`.

    Examples
    --------
    >>> ms = MultiSpeciesModel(lattice, species=2, q_per_species=3, steps=500)
    >>> ms.init_multispecies_dynamics()
    >>> ms.run()
    """

    dyn_UVclass = "multi_species"

    _c_program_name_template = "MultiSpeciesSimulator{}"
    _allowed_c_keys: tuple[str, ...] = ("C0",)

    def __init__(
        self,
        sg: "SignedGraph",
        species: int = 2,
        q_per_species: int | list[int] = 3,
        T: float = 1.0,
        interaction_matrix: NDArray | None = None,
        *,
        steps: int | None = None,
        simref: float | None = None,
        save_observables: bool = False,
        **kw: Any,
    ) -> None:
        dynpath = getattr(sg, 'path_data', None)
        if dynpath is not None:
            dynpath = Path(dynpath) / "multi_species"

        self.species = species

        if isinstance(q_per_species, int):
            self._q_per_species = [q_per_species] * species
        else:
            self._q_per_species = list(q_per_species)
            if len(self._q_per_species) != species:
                raise ValueError("q_per_species length must match species count.")

        q_total = max(self._q_per_species)  # used for VecDynSys.q

        super().__init__(
            sg,
            q=q_total,
            discrete=True,
            T=T,
            steps=steps,
            simref=simref,
            dynpath=dynpath,
            **kw,
        )

        if interaction_matrix is not None:
            self.interaction_matrix = np.asarray(interaction_matrix, dtype=np.float64)
            if self.interaction_matrix.shape != (species, species):
                raise ValueError(
                    f"Interaction matrix shape {self.interaction_matrix.shape} "
                    f"!= ({species}, {species})"
                )
        else:
            self.interaction_matrix = np.eye(species, dtype=np.float64)

        self.save_observables = save_observables
        self.ene: list[float] = []
        self.sini: NDArray | None = None

    @property
    def _state_shape(self) -> tuple[int, ...]:
        return (self.N, self.species)  # k components per node

    # ------------------------------------------------------------------
    # State initialisation override
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise multi-component state."""
        match self.ic:
            case "uniform" | "random" | "rand":
                self.s = np.zeros((self.N, self.species), dtype=np.int32)
                for k in range(self.species):
                    self.s[:, k] = np.random.randint(
                        0, self._q_per_species[k], size=self.N
                    )
            case "zero" | "zeros":
                self.s = np.zeros((self.N, self.species), dtype=np.int32)
            case "custom":
                if custom is None:
                    raise ValueError("Provide a custom state array.")
                self.s = np.asarray(custom, dtype=np.int32).copy()
                if self.s.shape != (self.N, self.species):
                    raise ValueError(
                        f"Shape mismatch: {self.s.shape} != ({self.N}, {self.species})"
                    )
            case _:
                raise ValueError(f"Unknown IC: {self.ic!r}")

    # ------------------------------------------------------------------
    # Energy
    # ------------------------------------------------------------------
    def energy(self, s: NDArray | None = None) -> float:
        """Multi-species energy."""
        if s is None:
            s = self.s
        A = self._get_adj_matrix()
        E = 0.0
        for i in range(self.N):
            for j in range(i + 1, self.N):
                if A[i, j] != 0:
                    for s1 in range(self.species):
                        for s2 in range(self.species):
                            E -= (
                                A[i, j]
                                * self.interaction_matrix[s1, s2]
                                * float(s[i, s1] == s[j, s2])
                            )
        return E

    def per_species_density(self, s: NDArray | None = None) -> list[NDArray]:
        """Fraction of nodes in each state, per species."""
        if s is None:
            s = self.s
        densities = []
        for k in range(self.species):
            counts = np.bincount(
                s[:, k].astype(int), minlength=self._q_per_species[k]
            )
            densities.append(counts / self.N)
        return densities

    # ------------------------------------------------------------------
    # Metropolis update
    # ------------------------------------------------------------------
    def ds1step(self, nd: int) -> None:
        """Single-node Metropolis: flip one random component."""
        A = self._get_adj_matrix()
        # Choose random species component to flip
        comp = np.random.randint(self.species)
        current = self.s[nd, comp]
        q_k = self._q_per_species[comp]

        # Propose new value != current
        proposal = np.random.randint(0, q_k - 1)
        if proposal >= current:
            proposal += 1

        # Energy change (only affected component matters)
        dE = 0.0
        for j in range(self.N):
            if A[nd, j] != 0:
                for s2 in range(self.species):
                    w = A[nd, j] * self.interaction_matrix[comp, s2]
                    dE += w * (
                        float(current == self.s[j, s2])
                        - float(proposal == self.s[j, s2])
                    )

        if dE <= 0 or (self.T > 0 and np.random.uniform() < np.exp(-dE / self.T)):
            self.s[nd, comp] = np.int32(proposal)

    # ------------------------------------------------------------------
    # Dynamics
    # ------------------------------------------------------------------
    def init_multispecies_dynamics(self, custom: Any = None, exName: str = "") -> None:
        self._check_c_backend_or_fallback()
        self.ene = []
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
            self.init_multispecies_dynamics()

    def run_py(self, verbose: bool = False, **kw: Any) -> None:
        """Monte Carlo sweeps."""
        for sweep in range(self.steps):
            nodes = np.random.permutation(self.N)
            for nd in nodes:
                self.ds1step(nd)
            if self.save_observables:
                self.ene.append(self.energy())
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
            f"{self.species}",
            f"{self.q}",
            f"{self.T:.12g}",
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
