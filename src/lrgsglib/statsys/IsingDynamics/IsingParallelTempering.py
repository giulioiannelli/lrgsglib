"""IsingParallelTempering — replica-exchange sampling (Layer-1 scheme).

K replicas of the same system evolve at the temperatures of a ladder, each
under the shared ``statsys._thermal.ThermalEngine`` (one compiled sweep per
rung — all closures read the live substrate state, so binding ``self.s`` to a
replica's array selects it with zero copying). Every ``steps_per_exchange``
sweeps, adjacent rungs attempt a configuration swap with the standard
replica-exchange Metropolis criterion ``min(1, e^{Δβ·ΔE})`` (M1 acceptance
applied to the pair; even–odd alternation so all adjacent pairs are tried).

Replica energies are tracked incrementally from the sweeps' ΔE — the swap
test needs no energy recompute. The standard observables record the COLDEST
rung (the physically interesting chain), once per exchange round; per-rung
traces live in :attr:`pt_energy` / :attr:`pt_magn` (intensive), and the
exchange diagnostic in :meth:`exchange_rate` (decision 25).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import tqdm

from .._thermal import THERMAL_RULES, ThermalEngine, resolve_tie_flip
from .defaults import (
    ISING_ORDER_DEFAULT,
    ISING_PT_LADDER_DEFAULT,
    ISING_PT_LADDERS,
    ISING_PT_N_EXCHANGES_DEFAULT,
    ISING_PT_N_REPLICAS_DEFAULT,
    ISING_PT_STEPS_PER_EXCHANGE_DEFAULT,
    ISING_PT_T_MAX_DEFAULT,
    ISING_PT_T_MIN_DEFAULT,
    ISING_RULE_DEFAULT,
    ISING_TIE_FLIP_DEFAULT,
    ISING_UPD_MODE_DEFAULT,
)
from .IsingBase import IsingBase

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["IsingParallelTempering", "generate_temperature_ladder"]


def generate_temperature_ladder(
    ladder: str,
    T_min: float,
    T_max: float,
    n_replicas: int,
) -> np.ndarray:
    """Build the replica-temperature ladder (ascending, rung 0 = coldest).

    'geometric' spaces the rungs uniformly in ``ln T`` (the standard choice
    for roughly-constant adjacent overlap); 'linear' uniformly in ``T``.
    """
    if ladder not in ISING_PT_LADDERS:
        raise ValueError(
            f"ladder must be one of {ISING_PT_LADDERS}; got {ladder!r}."
        )
    n = int(n_replicas)
    if n < 2:
        raise ValueError(f"n_replicas must be >= 2; got {n}.")
    if not 0.0 < T_min < T_max:
        raise ValueError(
            f"Need 0 < T_min < T_max; got T_min={T_min}, T_max={T_max}."
        )
    if ladder == "geometric":
        return T_min * (T_max / T_min) ** (np.arange(n) / (n - 1))
    return np.linspace(T_min, T_max, n)


class IsingParallelTempering(IsingBase):
    """Replica-exchange (parallel tempering) on the shared MCMC engine.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings.
    n_replicas : int
        Number of ladder rungs (K >= 2).
    T_min, T_max : float
        Ladder endpoints (``k_B = 1``); ``T_min`` must be > 0 (the swap
        criterion uses ``β = 1/T``).
    ladder : {'geometric', 'linear'}
        Rung spacing (see :func:`generate_temperature_ladder`).
    T_ladder : ndarray, optional
        Explicit ascending rung temperatures; overrides the named ladder.
    steps_per_exchange : int
        Monte Carlo sweeps per rung between exchange attempts.
    n_exchanges : int
        Number of exchange rounds; the recording horizon (``steps``) equals
        this, one record of the coldest rung per round.
    rule, order, tie_flip_p
        Same Layer-2 kinetics axes as :class:`IsingMetropolis`.
    **kwargs
        Forwarded to :class:`IsingBase`.

    After a run: :attr:`pt_energy` / :attr:`pt_magn` are ``(K, n_exchanges)``
    intensive per-rung traces, :attr:`pt_final_states` the final ``(K, N)``
    configurations (``self.s`` is the coldest), :attr:`T_ladder` the rungs,
    and :meth:`exchange_rate` the per-adjacent-pair swap acceptance.
    """

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        n_replicas: int = ISING_PT_N_REPLICAS_DEFAULT,
        T_min: float = ISING_PT_T_MIN_DEFAULT,
        T_max: float = ISING_PT_T_MAX_DEFAULT,
        ladder: str = ISING_PT_LADDER_DEFAULT,
        T_ladder: np.ndarray | None = None,
        steps_per_exchange: int = ISING_PT_STEPS_PER_EXCHANGE_DEFAULT,
        n_exchanges: int = ISING_PT_N_EXCHANGES_DEFAULT,
        rule: str = ISING_RULE_DEFAULT,
        order: str = ISING_ORDER_DEFAULT,
        tie_flip_p: float | str | None = None,
        **kwargs: Any,
    ) -> None:
        if rule not in THERMAL_RULES:
            raise ValueError(
                f"rule must be one of {THERMAL_RULES}; got {rule!r}."
            )
        if T_ladder is not None:
            rungs = np.asarray(T_ladder, dtype=np.float64)
            if rungs.ndim != 1 or len(rungs) < 2:
                raise ValueError(
                    "T_ladder must be a 1-D array with >= 2 rungs."
                )
            if np.any(rungs <= 0.0):
                raise ValueError(
                    "All ladder temperatures must be > 0 (the swap criterion "
                    f"uses beta = 1/T); got min={rungs.min()}."
                )
            if np.any(np.diff(rungs) <= 0.0):
                raise ValueError("T_ladder must be strictly ascending.")
        else:
            rungs = generate_temperature_ladder(
                ladder, T_min, T_max, n_replicas
            )

        self.T_ladder = rungs
        self.n_replicas = len(rungs)
        self.steps_per_exchange = int(steps_per_exchange)
        if self.steps_per_exchange < 1:
            raise ValueError(
                f"steps_per_exchange must be >= 1; got {steps_per_exchange}."
            )
        self.n_exchanges = int(n_exchanges)
        if self.n_exchanges < 1:
            raise ValueError(f"n_exchanges must be >= 1; got {n_exchanges}.")
        self.T = float(rungs[0])
        self.rule = rule
        self.upd_mode = ISING_UPD_MODE_DEFAULT
        self.order = order
        self.tie_flip_p = resolve_tie_flip(
            rule, tie_flip_p, default=ISING_TIE_FLIP_DEFAULT
        )
        self.pt_energy: np.ndarray | None = None
        self.pt_magn: np.ndarray | None = None
        self.pt_exchanges: list[tuple[int, int, int, bool]] = []
        self.pt_final_states: list[np.ndarray] | None = None

        super().__init__(sg, steps=self.n_exchanges, **kwargs)

        # Setup-time capability check on the coldest rung's engine.
        self._make_engine()

    def initialize_run_parameters(
        self, steps: int | None = None, simref: float | None = None
    ) -> None:
        if steps is not None or simref is not None:
            raise ValueError(
                "IsingParallelTempering derives its horizon from n_exchanges "
                "(one record per exchange round); configure n_exchanges/"
                "steps_per_exchange at construction instead of passing "
                "steps=/simref= to run()."
            )
        super().initialize_run_parameters()

    # ------------------------------------------------------------------
    # Native pybind backend: deliberately refused
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        raise NotImplementedError(
            "runlang='pb' is deliberately not wired for "
            "IsingParallelTempering: the native pt kernel's exchange "
            "criterion (LRGSG_pt.c::pt_attempt_exchange) is a tautology — "
            "every proposed swap is accepted (same bug as the legacy python "
            "path, and it feeds intensive energies into Δβ·ΔE). The python "
            "backend implements the correct min(1, e^{Δβ·ΔE}) criterion "
            "(ref M7); use runlang='py'."
        )

    # ------------------------------------------------------------------
    # Engine + sampling
    # ------------------------------------------------------------------
    def _make_engine(self, T: float | None = None) -> ThermalEngine:
        return ThermalEngine(
            self,
            rule=self.rule,
            upd_mode=self.upd_mode,
            T=float(self.T_ladder[0]) if T is None else float(T),
            order=self.order,
            tie_flip_p=self.tie_flip_p,
        )

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        n_rep = self.n_replicas
        n_sites = float(self.N)
        # Independent replica states; all engine closures read self.s, so
        # binding self.s = replicas[r] selects a replica without copying.
        replicas = [self.s.copy() for _ in range(n_rep)]
        sweeps = [
            self._make_engine(T=T_r).compile_step() for T_r in self.T_ladder
        ]
        energies = np.array(
            [self.total_energy(rep) for rep in replicas], dtype=np.float64
        )
        betas = 1.0 / self.T_ladder

        self.pt_energy = np.zeros((n_rep, self.n_exchanges))
        self.pt_magn = np.zeros((n_rep, self.n_exchanges))
        self.pt_exchanges = []

        rounds = (
            tqdm.tqdm(range(self.n_exchanges), desc=type(self).__name__)
            if tqdm_on
            else range(self.n_exchanges)
        )
        for ex_round in rounds:
            # Record the coldest rung at round start (standard observables).
            self.s = replicas[0]
            self.T = float(self.T_ladder[0])
            self._e_running = float(energies[0])
            self._record()

            # One block of sweeps per rung, energy tracked incrementally.
            for r in range(n_rep):
                self.s = replicas[r]
                for _ in range(self.steps_per_exchange):
                    energies[r] += sweeps[r]()
                self.pt_energy[r, ex_round] = energies[r] / n_sites
                self.pt_magn[r, ex_round] = float(
                    np.mean(replicas[r], dtype=np.float64)
                )

            # Even-odd alternating adjacent exchanges: swap the
            # configurations at rungs (i, i+1) with min(1, e^{Δβ·ΔE}),
            # Δβ = β_i − β_{i+1} > 0, ΔE = E_i − E_{i+1} — a lower-energy
            # configuration found on the hotter rung always moves down.
            start = ex_round % 2
            for i in range(start, n_rep - 1, 2):
                d_beta = betas[i] - betas[i + 1]
                d_E = energies[i] - energies[i + 1]
                arg = d_beta * d_E
                accepted = arg >= 0.0 or np.random.random() < np.exp(arg)
                if accepted:
                    replicas[i], replicas[i + 1] = (
                        replicas[i + 1],
                        replicas[i],
                    )
                    energies[i], energies[i + 1] = (
                        energies[i + 1],
                        energies[i],
                    )
                self.pt_exchanges.append((ex_round, i, i + 1, accepted))

        self.s = replicas[0]
        self.T = float(self.T_ladder[0])
        self._e_running = float(energies[0])
        self.pt_final_states = replicas

    def exchange_rate(self) -> np.ndarray:
        """Swap acceptance rate per adjacent rung pair (length K-1)."""
        n_pairs = self.n_replicas - 1
        rates = np.zeros(n_pairs)
        counts = np.zeros(n_pairs)
        for _, i, _, accepted in self.pt_exchanges:
            counts[i] += 1
            if accepted:
                rates[i] += 1
        return rates / np.maximum(counts, 1)
