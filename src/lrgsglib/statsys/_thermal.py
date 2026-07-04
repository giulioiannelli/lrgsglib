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
    "ClusterSubstrate",
    "ExchangeSubstrate",
    "resolve_acceptance",
    "resolve_tie_flip",
    "resolve_visit_order",
    "ThermalEngine",
    "WolffEngine",
    "SwendsenWangEngine",
    "KawasakiEngine",
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


@runtime_checkable
class ClusterSubstrate(ThermalSubstrate, Protocol):
    """Extension contract for cluster moves (Wolff / Swendsen–Wang).

    The model owns everything state- and Hamiltonian-specific: what a
    "satisfied" bond is (Ising: ``J_ij s_i s_j > 0``; Potts:
    ``J_ij δ(σ_i, σ_j) > 0``) and the energy change of flipping a cluster.
    The engines own the Fortuin–Kasteleyn activation probability
    ``p = 1 − e^{−2·sat/T}`` and the cluster growth/percolation itself.
    """

    def bond_satisfaction(self, nd: int) -> tuple[np.ndarray, np.ndarray]:
        """``(neighbor_indices, satisfactions)`` of site *nd*; a bond's
        satisfaction is its satisfied-coupling magnitude (> 0 activatable,
        <= 0 never activated)."""
        ...

    def edge_satisfactions(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(u, v, satisfactions)`` over every undirected edge once."""
        ...

    def flip_cluster(self, nodes: np.ndarray) -> float:
        """Flip the cluster *nodes* in place; return the energy change."""
        ...


@runtime_checkable
class ExchangeSubstrate(ThermalSubstrate, Protocol):
    """Extension contract for the pair-exchange (Kawasaki) move: dynamics
    conserving the order parameter by swapping two neighboring states."""

    def neighbor_indices(self, nd: int) -> np.ndarray:
        """Neighbor indices of site *nd*."""
        ...

    def swap_delta_E(self, u: int, v: int) -> float:
        """Energy change of swapping the states at *u* and *v*."""
        ...

    def commit_swap(self, u: int, v: int) -> None:
        """Swap the states at *u* and *v* (in place)."""
        ...


def resolve_tie_flip(
    rule: str,
    tie_flip_p: float | str | None,
    default: float | None = None,
) -> float | None:
    """Resolve the user-facing tie policy to a float (or None).

    Shared by every thermal scheme leaf (Metropolis / annealing / tempering):
    under ``rule='metropolis'`` the value may be a probability, a named preset
    from :data:`TIE_FLIP_PRESETS`, or ``None`` (→ *default*); the other rules
    bake their own ΔE=0 rate and reject an explicit value (plan §3b).
    """
    if rule == "metropolis":
        if tie_flip_p is None:
            return default
        if isinstance(tie_flip_p, str):
            try:
                return TIE_FLIP_PRESETS[tie_flip_p]
            except KeyError:
                raise ValueError(
                    f"Unknown tie_flip_p preset {tie_flip_p!r}; known "
                    f"presets: {tuple(TIE_FLIP_PRESETS)}."
                ) from None
        return float(tie_flip_p)
    if tie_flip_p is not None:
        raise ValueError(
            f"tie_flip_p is metropolis-only; rule {rule!r} bakes its own "
            "ΔE=0 rate."
        )
    return None


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


def resolve_visit_order(
    order: str, n_sites: int, rng: Any
) -> Callable[[], np.ndarray]:
    """Resolve an async site-visit ``order`` name (see :data:`ASYNC_ORDERS`)
    to a compiled zero-argument callable returning the sweep's site sequence.
    Shared by every single-site engine (thermal and exchange kinetics)."""
    if order not in ASYNC_ORDERS:
        raise ValueError(
            f"Unknown async order {order!r}; expected one of {ASYNC_ORDERS}."
        )
    if order == "random":

        def visit_order() -> np.ndarray:
            return rng.randint(0, n_sites, size=n_sites)

    elif order == "permutation":

        def visit_order() -> np.ndarray:
            sites = np.arange(n_sites)
            rng.shuffle(sites)
            return sites

    else:  # typewriter (legacy/debug)
        fixed = np.arange(n_sites)

        def visit_order() -> np.ndarray:
            return fixed

    return visit_order


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
        visit_order = resolve_visit_order(self.order, n_sites, rng)
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


def _cluster_p_add(T: float) -> Callable[[float], float]:
    """Compile the Fortuin–Kasteleyn bond-activation probability
    ``p(sat) = 1 − e^{−2·sat/T}`` for satisfied bonds (sat > 0); the T = 0
    limit activates every satisfied bond (refs M5/M6)."""
    if T > 0.0:
        inv_T = 1.0 / T

        def p_add(sat: float) -> float:
            return -math.expm1(-2.0 * sat * inv_T)

    else:

        def p_add(sat: float) -> float:
            return 1.0

    return p_add


class WolffEngine:
    """Wolff single-cluster kinetics (ref M6) over a :class:`ClusterSubstrate`.

    One "sweep" repeats single-cluster steps (grow from a random seed along
    satisfied bonds with the FK probability, flip the whole cluster — always
    accepted) until at least N sites have been flipped, so a sweep is
    comparable to N single-site attempts. Requires zero external field (the
    plain cluster rule breaks detailed balance under h; the ghost-spin
    extension is not wired) — the SCHEME validates that, since the field
    lives on the model, not here.

    Per-sweep flipped-site counts accumulate in :attr:`cluster_sizes`.
    """

    def __init__(
        self,
        substrate: ClusterSubstrate,
        *,
        T: float,
        rng: Any = None,
    ) -> None:
        if not isinstance(substrate, ClusterSubstrate):
            raise NotImplementedError(
                f"{type(substrate).__name__} does not implement the cluster "
                "substrate contract (bond_satisfaction / edge_satisfactions "
                "/ flip_cluster); cluster moves are unavailable."
            )
        if T < 0.0:
            raise ValueError(f"Temperature must be >= 0; got {T}.")
        self.substrate = substrate
        self.T = float(T)
        self.rng = rng if rng is not None else np.random
        self.cluster_sizes: list[int] = []

    def compile_step(self) -> Callable[[], float]:
        substrate = self.substrate
        rng = self.rng
        p_add = _cluster_p_add(self.T)
        n_sites = int(len(substrate.s))
        bond_satisfaction = substrate.bond_satisfaction
        flip_cluster = substrate.flip_cluster
        cluster_sizes = self.cluster_sizes
        in_cluster = np.zeros(n_sites, dtype=np.bool_)

        def cluster_step() -> tuple[float, int]:
            seed = int(rng.randint(0, n_sites))
            in_cluster.fill(False)
            in_cluster[seed] = True
            queue = [seed]
            head = 0
            while head < len(queue):
                node = queue[head]
                head += 1
                nbrs, sats = bond_satisfaction(node)
                for j, sat in zip(nbrs, sats):
                    j = int(j)
                    if in_cluster[j] or sat <= 0.0:
                        continue
                    if rng.random() < p_add(float(sat)):
                        in_cluster[j] = True
                        queue.append(j)
            nodes = np.asarray(queue, dtype=np.intp)
            return flip_cluster(nodes), len(queue)

        def sweep() -> float:
            dE_sweep = 0.0
            flipped = 0
            while flipped < n_sites:
                dE, size = cluster_step()
                dE_sweep += dE
                flipped += size
            cluster_sizes.append(flipped)
            return dE_sweep

        return sweep


class SwendsenWangEngine:
    """Swendsen–Wang multi-cluster kinetics (ref M5) over a
    :class:`ClusterSubstrate`.

    One sweep: activate every satisfied bond with the FK probability,
    partition the sites with union–find, flip each cluster with probability
    ½. Same zero-field requirement as :class:`WolffEngine` (validated by the
    scheme). Per-sweep cluster counts accumulate in :attr:`cluster_counts`.
    """

    def __init__(
        self,
        substrate: ClusterSubstrate,
        *,
        T: float,
        rng: Any = None,
    ) -> None:
        if not isinstance(substrate, ClusterSubstrate):
            raise NotImplementedError(
                f"{type(substrate).__name__} does not implement the cluster "
                "substrate contract (bond_satisfaction / edge_satisfactions "
                "/ flip_cluster); cluster moves are unavailable."
            )
        if T < 0.0:
            raise ValueError(f"Temperature must be >= 0; got {T}.")
        self.substrate = substrate
        self.T = float(T)
        self.rng = rng if rng is not None else np.random
        self.cluster_counts: list[int] = []

    def compile_step(self) -> Callable[[], float]:
        substrate = self.substrate
        rng = self.rng
        T = self.T
        n_sites = int(len(substrate.s))
        edge_satisfactions = substrate.edge_satisfactions
        flip_cluster = substrate.flip_cluster
        cluster_counts = self.cluster_counts

        def sweep() -> float:
            u_arr, v_arr, sat = edge_satisfactions()
            sat = np.asarray(sat, dtype=np.float64)
            # FK activation, vectorized: satisfied bonds only.
            if T > 0.0:
                p = -np.expm1(-2.0 * np.clip(sat, 0.0, None) / T)
            else:
                p = (sat > 0.0).astype(np.float64)
            active = rng.random(len(sat)) < p

            # Union-find with path halving + union by rank.
            parent = np.arange(n_sites)
            rank = np.zeros(n_sites, dtype=np.int32)

            def find(x: int) -> int:
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = int(parent[x])
                return x

            for u, v in zip(u_arr[active], v_arr[active]):
                ru, rv = find(int(u)), find(int(v))
                if ru == rv:
                    continue
                if rank[ru] < rank[rv]:
                    parent[ru] = rv
                elif rank[ru] > rank[rv]:
                    parent[rv] = ru
                else:
                    parent[rv] = ru
                    rank[ru] += 1

            roots = np.fromiter(
                (find(i) for i in range(n_sites)), dtype=np.intp, count=n_sites
            )
            unique_roots = np.unique(roots)
            flip_bit = rng.random(len(unique_roots)) < 0.5

            dE_sweep = 0.0
            for root, flip in zip(unique_roots, flip_bit):
                if flip:
                    nodes = np.flatnonzero(roots == root)
                    dE_sweep += flip_cluster(nodes)
            cluster_counts.append(len(unique_roots))
            return dE_sweep

        return sweep


class KawasakiEngine:
    """Pair-exchange (Kawasaki, ref M8) kinetics over an
    :class:`ExchangeSubstrate`: order-parameter-conserving dynamics.

    One sweep makes N exchange attempts: visit a site (same ``order`` axis
    as the single-site engine), pick a uniform random neighbor, and swap the
    two states with the chosen acceptance ``rule`` applied to the swap's ΔE.
    Attempts on equal states are null moves.
    """

    def __init__(
        self,
        substrate: ExchangeSubstrate,
        *,
        rule: str,
        T: float,
        order: str = "random",
        tie_flip_p: float | None = None,
        rng: Any = None,
    ) -> None:
        if not isinstance(substrate, ExchangeSubstrate):
            raise NotImplementedError(
                f"{type(substrate).__name__} does not implement the exchange "
                "substrate contract (neighbor_indices / swap_delta_E / "
                "commit_swap); the exchange move is unavailable."
            )
        if order not in ASYNC_ORDERS:
            raise ValueError(
                f"Unknown async order {order!r}; expected one of {ASYNC_ORDERS}."
            )
        # Validates rule / T / tie eagerly (setup-time capability error).
        resolve_acceptance(rule, T, tie_flip_p)
        self.substrate = substrate
        self.rule = rule
        self.T = float(T)
        self.order = order
        self.tie_flip_p = tie_flip_p
        self.rng = rng if rng is not None else np.random

    def compile_step(self) -> Callable[[], float]:
        substrate = self.substrate
        rng = self.rng
        p_accept = resolve_acceptance(self.rule, self.T, self.tie_flip_p)
        n_sites = int(len(substrate.s))
        visit_order = resolve_visit_order(self.order, n_sites, rng)
        neighbor_indices = substrate.neighbor_indices
        swap_delta_E = substrate.swap_delta_E
        commit_swap = substrate.commit_swap

        def sweep() -> float:
            s = substrate.s
            dE_sweep = 0.0
            for u in visit_order():
                u = int(u)
                nbrs = neighbor_indices(u)
                if len(nbrs) == 0:
                    continue
                v = int(nbrs[rng.randint(0, len(nbrs))])
                if s[u] == s[v]:
                    continue  # identity swap
                dE = swap_delta_E(u, v)
                p = p_accept(dE)
                if p >= 1.0 or rng.random() < p:
                    commit_swap(u, v)
                    dE_sweep += dE
            return dE_sweep

        return sweep
