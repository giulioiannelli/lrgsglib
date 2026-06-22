"""Incremental cluster-size-distribution tracker for the voter model.

A *cluster* (domain) is a connected component of the subgraph of **active**
edges. The activation predicate is ``b_ij * s_i * s_j > 0`` where the per-edge
sign ``b_ij`` selects the cluster mode (see ``defaults.py``):

* ``satisfied`` -- ``b_ij = sign(w_ij)``: domains of aligned signed opinion (a
  negative edge "agrees" when the spins are opposite). Reduces to same-opinion
  domains on an unsigned graph.
* ``rawspin``   -- ``b_ij = +1``: domains of equal raw spin value.

Because a single spin flip toggles the active state of **every** incident edge,
the partition is maintained without an O(N+E) recompute per step:

* **merges** (the flipped node now agrees with its formerly-frustrated
  neighbours) are union-find unions -- O(alpha) each, the "put two previously
  unspeaking clusters in contact" case;
* **splits** can only happen when the node had >= 2 agreeing neighbours that it
  now abandons. Only then does a *bounded* local scan run; if the abandoned
  neighbours re-connect within ``split_scan_cap`` nodes the domain is intact,
  otherwise exactly that one component is recomputed (rare; O(component)).

The brute-force ``cluster_components`` / ``cluster_size_distribution`` are the
correctness oracle the incremental ``ClusterTracker`` is validated against.

This module is the Python reference; the native C/pybind CTMC kernel mirrors the
same algorithm and is checked against this oracle for parity.
"""

from __future__ import annotations

from collections import Counter, deque

import numpy as np

from .defaults import CLUSTER_SPLIT_SCAN_CAP


def edge_sign_arrays(signs, cluster_mode: str):
    """Return the per-edge ``b_ij`` arrays for the requested ``cluster_mode``.

    ``signs`` is the ragged per-node array of ``sign(w_ij)`` (the convention used
    by ``VoterModel._gillespie_neighbors``). For ``satisfied`` the edge sign is
    ``sign(w_ij)``; for ``rawspin`` it is ``+1`` (edge sign ignored).
    """
    if cluster_mode == "satisfied":
        return [np.asarray(sg, dtype=np.int8) for sg in signs]
    if cluster_mode == "rawspin":
        return [np.ones(len(sg), dtype=np.int8) for sg in signs]
    raise ValueError(f"unknown cluster_mode={cluster_mode!r}")


# ---------------------------------------------------------------------------
# Brute-force oracle (full connected-components over active edges)
# ---------------------------------------------------------------------------
def cluster_components(s, idx, b) -> np.ndarray:
    """Component label per node via BFS over active edges (the oracle).

    ``idx[i]`` are node ``i``'s neighbours, ``b[i]`` the matching per-edge signs;
    edge ``(i, idx[i][k])`` is active iff ``b[i][k] * s[i] * s[idx[i][k]] > 0``.
    """
    n = len(idx)
    label = np.full(n, -1, dtype=np.int64)
    cur = 0
    for start in range(n):
        if label[start] != -1:
            continue
        label[start] = cur
        stack = [start]
        while stack:
            u = stack.pop()
            su = int(s[u])
            ju, bu = idx[u], b[u]
            for k in range(ju.size):
                v = int(ju[k])
                if label[v] == -1 and int(bu[k]) * su * int(s[v]) > 0:
                    label[v] = cur
                    stack.append(v)
        cur += 1
    return label


def cluster_size_distribution(s, idx, b) -> Counter:
    """Counter ``{cluster_size: number_of_clusters}`` (the oracle)."""
    label = cluster_components(s, idx, b)
    sizes = np.bincount(label)
    return Counter(int(x) for x in sizes)


# ---------------------------------------------------------------------------
# Incremental tracker
# ---------------------------------------------------------------------------
class ClusterTracker:
    """Maintain the active-edge partition under single-spin flips.

    Call :meth:`flip` with a node index **immediately before** the owning model
    flips ``s[i]`` (the tracker reads the pre-flip configuration to classify the
    incident edges). The configuration array ``s`` is shared by reference.
    """

    def __init__(self, s, idx, b, split_scan_cap: int = CLUSTER_SPLIT_SCAN_CAP):
        self.s = s
        self.idx = idx
        self.b = b
        self.n = len(idx)
        self.split_scan_cap = int(split_scan_cap)
        # Partition state: component id per node + explicit member sets, plus an
        # incrementally-maintained histogram {size: count_of_components}.
        self.comp_of = np.empty(self.n, dtype=np.int64)
        self.members: dict[int, set[int]] = {}
        self.size_count: Counter = Counter()
        self._next_id = 0
        self._rebuild_from_scratch()

    # -- public observables ------------------------------------------------
    def size_distribution(self) -> dict[int, int]:
        """Current ``{cluster_size: number_of_clusters}`` (positive counts)."""
        return {sz: c for sz, c in self.size_count.items() if c > 0}

    def sizes(self) -> list[int]:
        """Flat list of every current cluster size."""
        return [len(m) for m in self.members.values()]

    @property
    def n_clusters(self) -> int:
        return len(self.members)

    # -- construction ------------------------------------------------------
    def _new_comp(self, member_set: set[int]) -> int:
        cid = self._next_id
        self._next_id += 1
        self.members[cid] = member_set
        for u in member_set:
            self.comp_of[u] = cid
        self.size_count[len(member_set)] += 1
        return cid

    def _drop_comp(self, cid: int) -> None:
        self.size_count[len(self.members[cid])] -= 1
        del self.members[cid]

    def _rebuild_from_scratch(self) -> None:
        self.members.clear()
        self.size_count.clear()
        self._next_id = 0
        label = cluster_components(self.s, self.idx, self.b)
        groups: dict[int, set[int]] = {}
        for u, lab in enumerate(label):
            groups.setdefault(int(lab), set()).add(u)
        for member_set in groups.values():
            self._new_comp(member_set)

    # -- the incremental update -------------------------------------------
    def _active(self, i: int, k: int) -> bool:
        """Active state of edge ``(i, idx[i][k])`` under the *current* ``s``."""
        v = int(self.idx[i][k])
        return int(self.b[i][k]) * int(self.s[i]) * int(self.s[v]) > 0

    def flip(self, i: int) -> None:
        """Update the partition for a flip of node ``i`` (call BEFORE flipping).

        Pre-flip, ``A`` = neighbours sharing ``i``'s cluster (active edges, to be
        broken) and ``B`` = the rest (inactive edges, to become active). After the
        flip ``i`` leaves ``A`` (possible split of its old cluster) and joins the
        clusters of ``B`` (merges).
        """
        ji = self.idx[i]
        if ji.size == 0:
            return  # isolated node: always its own singleton, nothing changes
        A: list[int] = []
        B: list[int] = []
        si = int(self.s[i])
        for k in range(ji.size):
            v = int(ji[k])
            if int(self.b[i][k]) * si * int(self.s[v]) > 0:
                A.append(v)
            else:
                B.append(v)

        # --- 1. detach i from its old cluster, handling a possible split ---
        self._detach(i, A)

        # --- 2. merge i (now a fresh singleton) with the B-neighbours ------
        # NOTE: classification used pre-flip s; the actual flip of s[i] is done
        # by the caller right after this returns, so B is exactly the set of
        # neighbours i will agree with post-flip.
        if B:
            self._merge(i, B)

    def _detach(self, i: int, A: list[int]) -> None:
        """Remove ``i`` from its current cluster; resolve any resulting split."""
        cid = int(self.comp_of[i])
        comp = self.members[cid]
        if len(comp) == 1:
            # i was alone; it stays alone (will be re-homed by the merge step).
            return
        # i leaves -> at most a leaf/singleton removal cannot disconnect the rest.
        if len(A) <= 1:
            self._drop_comp(cid)
            comp.discard(i)
            self._new_comp(comp)              # remaining stays one cluster
            self._new_singleton(i)
            return
        # len(A) >= 2: a split is *possible*. Bounded scan from A[0].
        reached = self._bounded_scan(i, A)
        if reached is None:
            # No split detected within the cap (all A re-connected): the old
            # cluster simply loses i.
            self._drop_comp(cid)
            comp.discard(i)
            self._new_comp(comp)
            self._new_singleton(i)
            return
        # A genuine split (or cap exceeded): recompute exactly this component.
        self._drop_comp(cid)
        comp.discard(i)
        self._recompute_component(comp)
        self._new_singleton(i)

    def _new_singleton(self, i: int) -> None:
        self._new_comp({i})

    def _bounded_scan(self, i: int, A: list[int]):
        """BFS in the active subgraph (excluding ``i``) from ``A[0]``.

        Returns ``None`` if all of ``A`` are mutually re-connected (no split) --
        either reached within the cap, or the scan exhausted having met them all.
        Returns the (small) reached set if the scan exhausted *before* reaching
        every node of ``A`` (a confirmed split), or if the cap was hit (uncertain
        -> caller does an exact recompute).
        """
        targets = set(A)
        src = A[0]
        seen = {src}
        targets.discard(src)
        dq = deque([src])
        cap = self.split_scan_cap
        while dq:
            u = dq.popleft()
            ju, bu = self.idx[u], self.b[u]
            su = int(self.s[u])
            for k in range(ju.size):
                v = int(ju[k])
                if v == i or v in seen:
                    continue
                if int(bu[k]) * su * int(self.s[v]) > 0:
                    seen.add(v)
                    targets.discard(v)
                    if not targets:
                        return None       # all A re-connected -> no split
                    if len(seen) > cap:
                        return seen        # cap hit -> uncertain, force recompute
                    dq.append(v)
        # Scan exhausted: if some targets remain, A[0]'s piece is a real fragment.
        return seen if targets else None

    def _recompute_component(self, comp: set[int]) -> None:
        """Re-partition ``comp`` into active-edge components (exact, O(comp))."""
        seen: set[int] = set()
        for start in comp:
            if start in seen:
                continue
            frag = {start}
            seen.add(start)
            stack = [start]
            while stack:
                u = stack.pop()
                ju, bu = self.idx[u], self.b[u]
                su = int(self.s[u])
                for k in range(ju.size):
                    v = int(ju[k])
                    if v in comp and v not in seen \
                            and int(bu[k]) * su * int(self.s[v]) > 0:
                        seen.add(v)
                        frag.add(v)
                        stack.append(v)
            self._new_comp(frag)

    def _merge(self, i: int, B: list[int]) -> None:
        """Merge ``i``'s singleton with every distinct cluster in ``B``.

        Uses small-to-large: the largest participating cluster absorbs the rest.
        """
        roots = {int(self.comp_of[i])}
        roots.update(int(self.comp_of[v]) for v in B)
        if len(roots) == 1:
            return  # already one cluster (i and all of B coincide)
        keep = max(roots, key=lambda c: len(self.members[c]))
        target = self.members[keep]
        # detach keep's size from the histogram; re-add once after growth
        self.size_count[len(target)] -= 1
        for cid in roots:
            if cid == keep:
                continue
            src = self.members.pop(cid)
            self.size_count[len(src)] -= 1
            for u in src:
                self.comp_of[u] = keep
            target |= src
        self.size_count[len(target)] += 1
