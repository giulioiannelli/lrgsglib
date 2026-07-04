"""Generic finite-temperature MCMC engine (shared spin-model home).

Why this file exists (Ising plan §4, decision D-B2): the sweep / single-flip /
acceptance machinery is a *spin-model* layer, not a base-class layer. It is
written ONCE here against the substrate contract below and consumed by every
spin model — Ising (``BinDynSys``-based) today, Potts/XY/Heisenberg/
MultiSpecies (``VecDynSys``-based) as they modernize. Flow models (Kuramoto,
reaction-diffusion, coupled ODE) and the walker never import it.

Design contract (plan §5, configure-then-run): everything is resolved ONCE at
setup — the acceptance rule, the update order, the RNG — into plain callables;
the hot loop only calls them. No per-step branching on strings, ever.

Citations (``.agents/ref/00-REFERENCES.md`` §"MCMC / spin-dynamics kernels"):

- ``metropolis`` → M1 (Metropolis et al. 1953): accept with probability
  ``min(1, e^(−ΔE/T))`` (``k_B = 1`` units). The ``tie_flip_p`` knob at
  ``ΔE = 0`` is a project kinetic extension, not in M1.
- ``glauber`` → M2 (Glauber 1963): accept with probability
  ``1/(1 + e^(ΔE/T))``; its T→0 limit flips a tie with probability ½.
- ``heatbath`` → M4 (Newman & Barkema 1999): resample the site from its
  conditional Gibbs distribution — for a two-state (±1) substrate this
  coincides with the Glauber acceptance, so both names resolve to the same
  kernel here.
- ``upd_mode='gillespie'`` → M3 (Bortz–Kalos–Lebowitz 1975), Phase-2 wiring.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any, Protocol, runtime_checkable

import numpy as np

__all__ = [
    "THERMAL_RULES",
    "UPDATE_MODES",
    "ASYNC_ORDERS",
    "TIE_FLIP_PRESETS",
    "ThermalSubstrate",
    "resolve_acceptance",
    "ThermalEngine",
]

#: Acceptance rules (see the module-docstring citation map).
THERMAL_RULES: tuple[str, ...] = ("metropolis", "glauber", "heatbath")

#: Update schedules (Ising plan §9): random-sequential single-site,
#: checkerboard / sublattice-parallel, and rejection-free continuous-time
#: (BKL) kinetics. Phase 1 implements ``async``; ``sync``/``gillespie`` are
#: validated at setup and land with the native backends (Phase 2).
UPDATE_MODES: tuple[str, ...] = ("async", "sync", "gillespie")

#: Site-visit orders for ``upd_mode='async'`` (plan §3b): ``random`` picks a
#: site uniformly at random each update (with replacement — true asynchronous
#: kinetics, the default); ``permutation`` visits each site once per sweep in
#: a fresh random order; ``typewriter`` is the fixed 0..N−1 legacy/debug order.
ASYNC_ORDERS: tuple[str, ...] = ("random", "permutation", "typewriter")

#: Named presets for the Metropolis ΔE=0 tie probability (plan §3b):
#: ``standard`` reproduces M1 (e^0 = 1, ties always accepted); ``glauberT0``
#: matches M2's T→0 tie rate of ½; ``frozen`` makes ties absorbing (the
#: kinetic-Ising T=0 frozen rule).
TIE_FLIP_PRESETS: dict[str, float] = {
    "standard": 1.0,
    "glauberT0": 0.5,
    "frozen": 0.0,
}


@runtime_checkable
class ThermalSubstrate(Protocol):
    """The substrate contract (Ising plan §3) a spin model exposes.

    The engine never sees model internals — only these operations. This is
    what makes the engine write-once: ``IsingBase`` implements them on an int8
    ±1 state; a future ``PottsBase`` implements them on an integer-label
    state; the engine is unchanged.
    """

    #: Current state array (any dtype/shape the model chooses; ``len(s)`` is
    #: the number of sites).
    s: np.ndarray

    def propose_local(self, nd: int) -> Any:
        """Propose a new local configuration at site *nd* (e.g. the flipped
        spin, a new Potts label, a rotated XY angle)."""
        ...

    def delta_E(self, nd: int, proposal: Any) -> float:
        """Energy change if *proposal* were committed at site *nd* (includes
        the external-field term)."""
        ...

    def commit(self, nd: int, proposal: Any) -> None:
        """Apply *proposal* at site *nd* (in place)."""
        ...

    def order_parameter(self) -> float:
        """The model's scalar order parameter for the current state."""
        ...

    def bond_activation(self, u: int, v: int) -> float:
        """Bond-activation probability for cluster moves (Wolff / SW); models
        without cluster moves may raise ``NotImplementedError``."""
        ...


def _logistic_neg(x: float) -> float:
    """Overflow-safe ``1 / (1 + e^x)``."""
    if x >= 0.0:
        e = math.exp(-x)
        return e / (1.0 + e)
    return 1.0 / (1.0 + math.exp(x))


def resolve_acceptance(
    rule: str,
    T: float,
    tie_flip_p: float | None = None,
) -> Callable[[float], float]:
    """Resolve an acceptance rule to a compiled ``ΔE -> p_accept`` callable.

    Called ONCE at setup (never in the hot loop): the rule name, the
    temperature regime (T > 0 vs the T = 0 limit) and the tie policy are all
    branched on here, so the returned callable is branch-minimal.

    Parameters
    ----------
    rule : {'metropolis', 'glauber', 'heatbath'}
        Acceptance rule (citations in the module docstring).
    T : float
        Temperature in ``k_B = 1`` units; ``T = 0`` selects the exact
        zero-temperature limit of the rule.
    tie_flip_p : float, optional
        Metropolis-only ΔE=0 acceptance probability (see
        :data:`TIE_FLIP_PRESETS`). Omitted → 1.0, the standard M1 behaviour.
        ``glauber``/``heatbath`` bake their own tie rate (½) and reject this
        parameter (plan §3b).
    """
    if rule not in THERMAL_RULES:
        raise ValueError(
            f"Unknown acceptance rule {rule!r}; expected one of {THERMAL_RULES}."
        )
    T = float(T)
    if T < 0.0:
        raise ValueError(f"Temperature must be >= 0; got {T}.")

    if rule == "metropolis":
        tie = (
            TIE_FLIP_PRESETS["standard"]
            if tie_flip_p is None
            else float(tie_flip_p)
        )
        if not 0.0 <= tie <= 1.0:
            raise ValueError(f"tie_flip_p must be in [0, 1]; got {tie}.")
        if T > 0.0:

            def p_accept(dE: float) -> float:
                if dE < 0.0:
                    return 1.0
                if dE == 0.0:
                    return tie
                return math.exp(-dE / T)

        else:

            def p_accept(dE: float) -> float:
                if dE < 0.0:
                    return 1.0
                if dE == 0.0:
                    return tie
                return 0.0

        return p_accept

    # glauber / heatbath (identical kernel on a two-state substrate; M2/M4)
    if tie_flip_p is not None:
        raise ValueError(
            f"tie_flip_p is a metropolis-only knob; rule {rule!r} bakes its "
            "own ΔE=0 rate (1/2)."
        )
    if T > 0.0:

        def p_accept(dE: float) -> float:
            return _logistic_neg(dE / T)

    else:

        def p_accept(dE: float) -> float:
            if dE < 0.0:
                return 1.0
            if dE == 0.0:
                return 0.5
            return 0.0

    return p_accept


class ThermalEngine:
    """Configure-then-run MCMC driver over a :class:`ThermalSubstrate`.

    Everything is resolved to callables in :meth:`compile_step`; the returned
    sweep function performs ``N`` single-site update attempts with zero
    per-step branching on configuration. Scheme classes (``IsingMetropolis``,
    ``IsingSimulatedAnnealing``, ...) own the temperature *policy* (fixed T,
    T(t) schedule, replica ladder) and drive this engine; the engine owns the
    *mechanics* of a sweep.

    Parameters
    ----------
    substrate : ThermalSubstrate
        The spin model (state + propose/ΔE/commit).
    rule : {'metropolis', 'glauber', 'heatbath'}
        Acceptance rule.
    upd_mode : {'async', 'sync', 'gillespie'}
        Update schedule. Phase 1 implements ``async``; the others raise
        ``NotImplementedError`` here at setup (hard capability error, no
        silent fallback).
    T : float
        Temperature (``k_B = 1``); mutable between sweeps via :attr:`T` +
        :meth:`compile_step` re-resolution (annealing recompiles per stage).
    order : {'random', 'permutation', 'typewriter'}
        Site-visit order for ``async`` (default ``'random'``, plan §3b).
    tie_flip_p : float, optional
        Metropolis-only tie probability (see :func:`resolve_acceptance`).
    rng : numpy RandomState-like, optional
        Random source with ``randint``/``shuffle``/``random`` — defaults to
        the global ``numpy.random`` module, which ``DynSys`` seeds, so seeded
        runs reproduce.
    """

    def __init__(
        self,
        substrate: ThermalSubstrate,
        *,
        rule: str,
        upd_mode: str,
        T: float,
        order: str = "random",
        tie_flip_p: float | None = None,
        rng: Any = None,
    ) -> None:
        if upd_mode not in UPDATE_MODES:
            raise ValueError(
                f"Unknown upd_mode {upd_mode!r}; expected one of {UPDATE_MODES}."
            )
        if upd_mode != "async":
            raise NotImplementedError(
                f"upd_mode={upd_mode!r} is not wired yet (Phase 2: 'sync' with "
                "the native backends, 'gillespie' via the shared CTMC core, "
                "ref M3). Use upd_mode='async'."
            )
        if order not in ASYNC_ORDERS:
            raise ValueError(
                f"Unknown async order {order!r}; expected one of {ASYNC_ORDERS}."
            )
        # Validates rule / T / tie_flip_p eagerly (setup-time capability error).
        resolve_acceptance(rule, T, tie_flip_p)
        self.substrate = substrate
        self.rule = rule
        self.upd_mode = upd_mode
        self.T = float(T)
        self.order = order
        self.tie_flip_p = tie_flip_p
        self.rng = rng if rng is not None else np.random

    def compile_step(self) -> Callable[[], float]:
        """Resolve (rule × T × order) ONCE into a sweep callable.

        The returned zero-argument function performs ``N`` single-site update
        attempts and returns the total energy change of the sweep, so the
        caller can track the running energy incrementally instead of
        recomputing the full Hamiltonian every sweep.
        """
        substrate = self.substrate
        rng = self.rng
        p_accept = resolve_acceptance(self.rule, self.T, self.tie_flip_p)
        n_sites = int(len(substrate.s))

        if self.order == "random":

            def visit_order() -> np.ndarray:
                return rng.randint(0, n_sites, size=n_sites)

        elif self.order == "permutation":

            def visit_order() -> np.ndarray:
                sites = np.arange(n_sites)
                rng.shuffle(sites)
                return sites

        else:  # typewriter (legacy/debug)
            fixed = np.arange(n_sites)

            def visit_order() -> np.ndarray:
                return fixed

        propose = substrate.propose_local
        delta_E = substrate.delta_E
        commit = substrate.commit

        def sweep() -> float:
            dE_sweep = 0.0
            for nd in visit_order():
                nd = int(nd)
                proposal = propose(nd)
                dE = delta_E(nd, proposal)
                p = p_accept(dE)
                if p >= 1.0 or rng.random() < p:
                    commit(nd, proposal)
                    dE_sweep += dE
            return dE_sweep

        return sweep
