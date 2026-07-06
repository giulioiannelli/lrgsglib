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
- ``upd_mode='gillespie'`` → M3 (Bortz–Kalos–Lebowitz 1975): rejection-free
  continuous-time kinetics from the per-site rate table.
- ``upd_mode='sync'`` → M4 (Newman & Barkema 1999): sublattice-parallel
  updates — color the interaction graph, update one independent set at a
  time from the frozen state (checkerboard on a bipartite lattice).
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

#: Update schedules (Ising plan §9): ``async`` = random-sequential
#: single-site attempts; ``sync`` = sublattice-parallel (greedy-color the
#: interaction graph once, then update one independent set at a time from the
#: frozen state — ref M4); ``gillespie`` = rejection-free continuous-time BKL
#: kinetics (ref M3; one recorded step = one unit of time = one
#: sweep-equivalent). ``sync``/``gillespie`` need ``neighbor_indices`` on the
#: substrate and take no async ``order``.
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
        Update schedule (see :data:`UPDATE_MODES`): random-sequential
        attempts, sublattice-parallel synchronous updates (ref M4), or
        rejection-free continuous-time BKL kinetics (ref M3). ``sync`` and
        ``gillespie`` require ``substrate.neighbor_indices`` and take no
        async ``order`` (setup-time capability errors, no silent fallback).
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
            # sync (sublattice coloring) and gillespie (BKL rate refresh)
            # both walk the adjacency, so the substrate must expose it. The
            # arrays behind it are only touched at compile_step time (after
            # the model's init), so hasattr is the right setup-time check.
            if not hasattr(substrate, "neighbor_indices"):
                raise NotImplementedError(
                    f"upd_mode={upd_mode!r} needs the substrate to expose "
                    f"neighbor_indices(nd); {type(substrate).__name__} "
                    "does not."
                )
            if order != "random":
                raise ValueError(
                    "order applies to upd_mode='async' only; "
                    f"upd_mode={upd_mode!r} has no site-visit order."
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
        #: Greedy color classes (sync mode only), cached across recompiles.
        self.color_classes: list[np.ndarray] | None = None

    def compile_step(self) -> Callable[[], float]:
        """Resolve (rule × T × upd_mode × order) ONCE into a sweep callable.

        The returned zero-argument function advances the system by one
        sweep-equivalent (N update attempts for ``async``/``sync``, one unit
        of continuous time for ``gillespie``) and returns the total energy
        change, so the caller can track the running energy incrementally
        instead of recomputing the full Hamiltonian every sweep.
        """
        if self.upd_mode == "sync":
            return self._compile_sync()
        if self.upd_mode == "gillespie":
            return self._compile_gillespie()
        return self._compile_async()

    def _compile_async(self) -> Callable[[], float]:
        """Random-sequential kinetics: N single-site attempts per sweep."""
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

    def _greedy_color_classes(self) -> list[np.ndarray]:
        """Greedy-color the interaction graph once: no two sites in a class
        are adjacent, so their ΔE's from the frozen state are simultaneously
        valid and the class updates as an independent block (ref M4). On a
        bipartite lattice this recovers the checkerboard (2 classes)."""
        substrate = self.substrate
        n_sites = int(len(substrate.s))
        # Presence was validated at setup; the base protocol stays minimal.
        neighbor_indices = getattr(substrate, "neighbor_indices")
        colors = np.full(n_sites, -1, dtype=np.int64)
        for i in range(n_sites):
            nbrs = np.asarray(neighbor_indices(i), dtype=np.intp)
            used = {int(c) for c in colors[nbrs] if c >= 0}
            c = 0
            while c in used:
                c += 1
            colors[i] = c
        return [
            np.flatnonzero(colors == c) for c in range(int(colors.max()) + 1)
        ]

    def _compile_sync(self) -> Callable[[], float]:
        """Sublattice-parallel synchronous kinetics (ref M4): per sweep,
        visit each color class in turn — compute every ΔE in the class from
        the frozen state first, then accept/commit independently. Each block
        preserves the Gibbs measure, so the composition does too."""
        substrate = self.substrate
        rng = self.rng
        p_accept = resolve_acceptance(self.rule, self.T, self.tie_flip_p)
        if self.color_classes is None:
            self.color_classes = self._greedy_color_classes()
        classes = [[int(nd) for nd in cls] for cls in self.color_classes]
        propose = substrate.propose_local
        delta_E = substrate.delta_E
        commit = substrate.commit

        def sweep() -> float:
            dE_sweep = 0.0
            for cls in classes:
                # Pass 1: proposals + ΔE for the whole class, state frozen.
                proposals = [propose(nd) for nd in cls]
                dEs = [delta_E(nd, prop) for nd, prop in zip(cls, proposals)]
                # Pass 2: independent accept/commit (no intra-class bonds,
                # so the pass-1 ΔE's stay exact as commits land).
                for nd, prop, dE in zip(cls, proposals, dEs):
                    p = p_accept(dE)
                    if p >= 1.0 or rng.random() < p:
                        commit(nd, prop)
                        dE_sweep += dE
            return dE_sweep

        return sweep

    def _compile_gillespie(self) -> Callable[[], float]:
        """Rejection-free continuous-time (BKL, ref M3) kinetics.

        Each site attempts at unit rate, so site *i* fires at
        ``r_i = p_accept(ΔE_i)`` and one unit of simulated time is one
        sweep-equivalent (N attempts on average). A sweep call advances time
        by exactly 1.0: draw waiting times ``Δt ~ Exp(Σ r)``, pick the event
        site ∝ ``r_i``, commit it (always accepted), refresh the rates of the
        site and its neighbors; the exponential's memorylessness lets the
        step truncate at the record boundary without executing the overshoot.
        ``Σ r = 0`` is an absorbing state — the remaining time just elapses.
        Assumes the local proposal is deterministic (two-state substrate);
        multi-state models need per-channel rate tables.
        """
        substrate = self.substrate
        rng = self.rng
        p_accept = resolve_acceptance(self.rule, self.T, self.tie_flip_p)
        n_sites = int(len(substrate.s))
        propose = substrate.propose_local
        delta_E = substrate.delta_E
        commit = substrate.commit
        # Presence was validated at setup; the base protocol stays minimal.
        neighbor_indices = getattr(substrate, "neighbor_indices")
        rates = np.zeros(n_sites, dtype=np.float64)
        dE_all = np.zeros(n_sites, dtype=np.float64)

        def refresh(nd: int) -> None:
            dE = delta_E(nd, propose(nd))
            dE_all[nd] = dE
            rates[nd] = p_accept(dE)

        def sweep() -> float:
            # Rebuild the rate table at the step boundary (O(N), the same
            # order as recording an observable): robust to any state change
            # between recorded steps (re-init, annealing recompile).
            for nd in range(n_sites):
                refresh(nd)
            dE_sweep = 0.0
            t = 0.0
            while True:
                R = float(rates.sum())
                if R <= 0.0:
                    break  # absorbing: no event can fire
                t += -math.log(1.0 - rng.random()) / R
                if t >= 1.0:
                    break  # overshoot past the record boundary (memoryless)
                cum = np.cumsum(rates)
                k = int(np.searchsorted(cum, rng.random() * R, side="right"))
                if k >= n_sites:  # float-roundoff guard on the last bin
                    k = n_sites - 1
                commit(k, propose(k))
                dE_sweep += dE_all[k]
                refresh(k)
                for j in neighbor_indices(k):
                    refresh(int(j))
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
