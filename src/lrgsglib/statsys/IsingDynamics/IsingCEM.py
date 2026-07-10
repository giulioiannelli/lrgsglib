"""IsingCEM — cross-entropy ground-state optimizer (Layer-1 scheme).

The fourth irreducible algorithm of the Ising family (decision 4): a
population optimizer with NO temperature — which is why ``T`` lives on the
thermal schemes and not on :class:`IsingBase`. It searches the spectral
coefficient space: candidate configurations are ``s = sign(c · V)`` with
``V`` the ``(M, N)`` basis of the lowest signed-Laplacian eigenvectors
(``statsys._spectral``), and the coefficient distribution ``N(μ, σ)`` is
iteratively refit to the elite fraction of each population (cross-entropy
update with exponential smoothing), optionally polishing candidates with a
zero-temperature greedy quench on the model's own CSR couplings.

Unlike the legacy ``cem_spectral_sampling``, the objective here is the FULL
§3b Hamiltonian ``H = −Σ J_ij s_i s_j − Σ h_i s_i`` via ``total_energy`` —
coupling_norm and the external field are honored, and the greedy polish
minimizes the same H it reports (plan §3b consistency by construction).

Observable semantics: one record per CEM iteration, tracing the INCUMBENT
(the restart's best-so-far configuration) — energy intensive, as everywhere.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import tqdm

from .._spectral import eigvec_subspace_matrix, greedy_polish_csr
from .defaults import (
    ISING_CEM_ELITE_FRAC_DEFAULT,
    ISING_CEM_GREEDY_DEFAULT,
    ISING_CEM_GREEDY_SWEEPS_DEFAULT,
    ISING_CEM_INIT_SIGMA_DEFAULT,
    ISING_CEM_ITER_DEFAULT,
    ISING_CEM_POP_SIZE_DEFAULT,
    ISING_CEM_RESTARTS_DEFAULT,
    ISING_CEM_SIGMA_CEILING_DEFAULT,
    ISING_CEM_SIGMA_FLOOR_DEFAULT,
    ISING_CEM_SMOOTHING_DEFAULT,
    ISING_SPECTRAL_N_MODES_DEFAULT,
    ISING_SPECTRAL_POLISH_DEFAULT,
    ISING_SPECTRAL_POLISH_SWEEPS_DEFAULT,
)
from .IsingBase import IsingBase

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["IsingCEM"]

# Additive floor keeping the refit σ strictly positive when an elite set
# happens to be degenerate (all-equal coefficients).
_SIGMA_EPS: float = 1e-12


class IsingCEM(IsingBase):
    """Cross-entropy method in the eigenvector-coefficient space.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings.
    n_modes : int
        Subspace dimension M (capped at N).
    pop_size : int
        Candidates per CEM iteration (K).
    elite_frac : float
        Fraction of the population refitting the distribution.
    n_iter : int
        CEM iterations per restart.
    restarts : int
        Independent restarts; the global best is kept.
    init_sigma : float
        Initial per-mode σ of the coefficient Gaussian.
    smoothing : float
        Exponential smoothing α of the (μ, σ) update:
        ``new = α·old + (1−α)·refit``.
    sigma_floor, sigma_ceiling : float
        Clip bounds on σ after each refit.
    greedy, greedy_sweeps
        Zero-temperature greedy polish of every candidate.
    polish, polish_sweeps
        Final greedy polish of each restart's best.
    **kwargs
        Forwarded to :class:`IsingBase`.

    After a run: :attr:`cem_best_spins` / :attr:`cem_best_energy` (intensive)
    / :attr:`cem_best_coeffs` hold the global best (``self.s`` is set to it),
    :attr:`cem_restart_energies` the per-restart bests (intensive).
    """

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        n_modes: int = ISING_SPECTRAL_N_MODES_DEFAULT,
        pop_size: int = ISING_CEM_POP_SIZE_DEFAULT,
        elite_frac: float = ISING_CEM_ELITE_FRAC_DEFAULT,
        n_iter: int = ISING_CEM_ITER_DEFAULT,
        restarts: int = ISING_CEM_RESTARTS_DEFAULT,
        init_sigma: float = ISING_CEM_INIT_SIGMA_DEFAULT,
        smoothing: float = ISING_CEM_SMOOTHING_DEFAULT,
        sigma_floor: float = ISING_CEM_SIGMA_FLOOR_DEFAULT,
        sigma_ceiling: float = ISING_CEM_SIGMA_CEILING_DEFAULT,
        greedy: bool = ISING_CEM_GREEDY_DEFAULT,
        greedy_sweeps: int = ISING_CEM_GREEDY_SWEEPS_DEFAULT,
        polish: bool = ISING_SPECTRAL_POLISH_DEFAULT,
        polish_sweeps: int = ISING_SPECTRAL_POLISH_SWEEPS_DEFAULT,
        **kwargs: Any,
    ) -> None:
        self.n_modes = int(n_modes)
        self.pop_size = int(pop_size)
        self.elite_frac = float(elite_frac)
        self.n_iter = int(n_iter)
        self.restarts = int(restarts)
        self.init_sigma = float(init_sigma)
        self.smoothing = float(smoothing)
        self.sigma_floor = float(sigma_floor)
        self.sigma_ceiling = float(sigma_ceiling)
        self.greedy = bool(greedy)
        self.greedy_sweeps = int(greedy_sweeps)
        self.polish = bool(polish)
        self.polish_sweeps = int(polish_sweeps)
        for name in ("pop_size", "n_iter", "restarts", "n_modes"):
            if getattr(self, name) < 1:
                raise ValueError(
                    f"{name} must be >= 1; got {getattr(self, name)}."
                )
        if not 0.0 < self.elite_frac <= 1.0:
            raise ValueError(
                f"elite_frac must be in (0, 1]; got {self.elite_frac}."
            )
        if not 0.0 <= self.smoothing < 1.0:
            raise ValueError(
                f"smoothing must be in [0, 1); got {self.smoothing}."
            )

        self.cem_best_spins: np.ndarray | None = None
        self.cem_best_energy: float | None = None
        self.cem_best_coeffs: np.ndarray | None = None
        self.cem_restart_energies: np.ndarray | None = None

        # One record per CEM iteration (the incumbent trace).
        super().__init__(sg, steps=self.restarts * self.n_iter, **kwargs)

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C)
    # ------------------------------------------------------------------
    _sch_token = "cem"

    def _ns_token_value(self) -> int | None:
        return None  # horizon derived from the scheme numerics

    def _physics_tokens(self) -> list:
        return [
            ("nrst", int(self.restarts)),
            ("nit", int(self.n_iter)),
            ("sig", float(self.init_sigma)),
        ]

    def initialize_run_parameters(
        self, steps: int | None = None, simref: float | None = None
    ) -> None:
        if steps is not None or simref is not None:
            raise ValueError(
                "IsingCEM derives its horizon from restarts x n_iter (one "
                "record per CEM iteration); configure those at construction "
                "instead of passing steps=/simref= to run()."
            )
        super().initialize_run_parameters()

    # ------------------------------------------------------------------
    # Sampling (no ThermalEngine — population optimizer)
    # ------------------------------------------------------------------
    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        n_modes = min(self.n_modes, self.N)
        V = eigvec_subspace_matrix(self.sg, n_modes)  # (M, N)
        self._ensure_couplings()
        nbr_idx, nbr_J, nbr_ptr = self._nbr_idx, self._nbr_J, self._nbr_ptr
        field = np.asarray(self.field, dtype=np.float64)

        K = self.pop_size
        n_elite = max(1, int(self.elite_frac * K))

        global_best_E = np.inf
        global_best_spins: np.ndarray | None = None
        global_best_coeffs: np.ndarray | None = None
        restart_E = np.zeros(self.restarts, dtype=np.float64)

        outer = (
            tqdm.tqdm(range(self.restarts), desc=type(self).__name__)
            if tqdm_on
            else range(self.restarts)
        )
        for r in outer:
            mu = np.zeros(n_modes, dtype=np.float64)
            sigma = np.full(n_modes, self.init_sigma, dtype=np.float64)
            best_E = np.inf
            best_spins: np.ndarray | None = None
            best_coeffs: np.ndarray | None = None

            for _ in range(self.n_iter):
                # 1. Sample K coefficient vectors; map to ±1 candidates.
                coeffs_pop = np.random.normal(
                    loc=mu, scale=sigma, size=(K, n_modes)
                )
                spins_pop = np.sign(coeffs_pop @ V).astype(np.int8)
                spins_pop[spins_pop == 0] = 1

                # 2. Optional zero-temperature greedy quench per candidate.
                if self.greedy:
                    for k in range(K):
                        spins_pop[k] = greedy_polish_csr(
                            spins_pop[k],
                            nbr_idx,
                            nbr_J,
                            nbr_ptr,
                            field=field,
                            max_sweeps=self.greedy_sweeps,
                        )

                # 3. Evaluate the full §3b Hamiltonian; select elites.
                energies = np.array(
                    [self.total_energy(spins_pop[k]) for k in range(K)]
                )
                elite_idx = np.argsort(energies)[:n_elite]

                # 4. Refit N(mu, sigma) to the elites, smoothed and clipped.
                elite = coeffs_pop[elite_idx]
                alpha = self.smoothing
                mu = alpha * mu + (1.0 - alpha) * elite.mean(axis=0)
                sigma = alpha * sigma + (1.0 - alpha) * (
                    elite.std(axis=0) + _SIGMA_EPS
                )
                sigma = np.clip(sigma, self.sigma_floor, self.sigma_ceiling)

                best_k = int(elite_idx[0])
                if energies[best_k] < best_E:
                    best_E = float(energies[best_k])
                    best_spins = spins_pop[best_k].copy()
                    best_coeffs = coeffs_pop[best_k].copy()

                # Incumbent trace: one record per CEM iteration.
                self.s = best_spins
                self._e_running = best_E
                self._record()

            # Optional final polish of the restart best.
            if self.polish and best_spins is not None:
                polished = greedy_polish_csr(
                    best_spins,
                    nbr_idx,
                    nbr_J,
                    nbr_ptr,
                    field=field,
                    max_sweeps=self.polish_sweeps,
                )
                E_pol = self.total_energy(polished)
                if E_pol < best_E:
                    best_E, best_spins = float(E_pol), polished

            restart_E[r] = best_E
            if best_E < global_best_E:
                global_best_E = best_E
                global_best_spins = best_spins
                global_best_coeffs = best_coeffs

        if global_best_spins is not None:
            self.s = global_best_spins.astype(np.int8)
        self._e_running = global_best_E
        self.cem_best_spins = global_best_spins
        self.cem_best_energy = global_best_E / float(self.N)
        self.cem_best_coeffs = global_best_coeffs
        self.cem_restart_energies = restart_E / float(self.N)
