"""Pure-numpy walker kernels for the SignedWalker family.

Engine-agnostic: operates on dense ``(N, z)`` neighbour / sign tables
extracted from any :class:`SignedGraphProtocol` instance via
:func:`signed_lattice_tables`.

Three rules, each in its own kernel:

* :func:`_absorb_kernel` — walker dies when it moves across an edge
  whose "kill-mask" entry is True.
* :func:`_kill_kernel`   — walker dies with probability ``1/n_neg`` at
  any node that has at least one negative incident edge.
* :func:`_sticky_kernel` — walker's move probability halves each time
  it crosses a frustrated edge; dies once it has crossed
  ``max_n_cross`` of them. Real-time is advanced by a geometric
  waiting time so idle steps are skipped.

X-node behaviour is captured by the ``kill_mask`` argument:

* ``'reflect'`` mode (library default) — only **frustrated** crossings
  kill; **unfrustrated** negative edges (incident to an X-node) act as
  reflective barriers (walker stays at its current site).
* ``'absorb'`` mode (notebook-default) — **any** negative crossing
  kills; unfrustrated negative edges are also absorbing.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from scipy.sparse import issparse


def signed_lattice_tables(sg: Any):
    """Dense ``(N, z)`` neighbour / sign tables for a z-regular signed graph.

    Extracts from ``sg.adj`` (scipy sparse for NX, dense ndarray for GT)
    an ``(N, z)`` integer neighbour index table, an ``(N, z)`` boolean
    negative-edge table, and an ``(N, z)`` boolean frustrated-edge table.

    A node is an **X-error node** iff all its incident edges are
    negative. A negative edge is **frustrated** iff neither endpoint is
    an X-error node.

    Parameters
    ----------
    sg : SignedGraphProtocol
        Must expose ``.adj`` (scipy-sparse or ndarray with ±1 entries)
        and ``.N`` (node count).

    Returns
    -------
    nbrs : (N, z) int64
    neg  : (N, z) bool
    frust: (N, z) bool
    """
    A = sg.adj
    if issparse(A):
        A = A.tocsr()
        N = A.shape[0]
        deg = np.diff(A.indptr)
        z = int(deg[0])
        if not np.all(deg == z):
            raise ValueError(
                f"signed_lattice_tables requires a z-regular graph; "
                f"got degrees in [{deg.min()}, {deg.max()}]."
            )
        nbrs = A.indices.reshape(N, z).astype(np.int64)
        signs = A.data.reshape(N, z).astype(np.int8)
    else:
        A_dense = np.asarray(A)
        N = A_dense.shape[0]
        mask = A_dense != 0
        deg = mask.sum(axis=1)
        z = int(deg[0])
        if not np.all(deg == z):
            raise ValueError(
                f"signed_lattice_tables requires a z-regular graph; "
                f"got degrees in [{deg.min()}, {deg.max()}]."
            )
        nbrs = np.empty((N, z), dtype=np.int64)
        signs = np.empty((N, z), dtype=np.int8)
        for i in range(N):
            js = np.flatnonzero(mask[i])
            nbrs[i] = js
            signs[i] = np.sign(A_dense[i, js]).astype(np.int8)
    neg = signs == -1
    antiferro = neg.all(axis=1)
    frust = neg & ~antiferro[:, None] & ~antiferro[nbrs]
    return nbrs, neg, frust


def _kill_masks(
    neg: np.ndarray,
    frust: np.ndarray,
    x_node_behavior: str,
):
    """Return ``(kill_mask, reflect_mask)`` tables for a given X-node mode.

    * ``x_node_behavior='reflect'`` (library default)
      - kill when moving across a frustrated edge
      - reflect (stay put) when moving across a non-frustrated negative
        edge (unfrustrated neighbours of X-nodes)
    * ``x_node_behavior='absorb'`` (notebook-default)
      - kill when moving across ANY negative edge
      - never reflect
    """
    if x_node_behavior == 'reflect':
        kill_mask = frust
        reflect_mask = neg & ~frust
    elif x_node_behavior == 'absorb':
        kill_mask = neg
        reflect_mask = np.zeros_like(neg)
    else:
        raise ValueError(
            f"x_node_behavior must be 'reflect' or 'absorb', "
            f"got {x_node_behavior!r}."
        )
    return kill_mask, reflect_mask


def _init_walker_state(
    n_walkers: int,
    pos_idx: np.ndarray,
    N: int,
    stop_thresh: int,
    store_trajectory: bool,
    store_per_walker_visits: bool,
):
    """Allocate and pre-fill the per-walker bookkeeping arrays."""
    wid = np.arange(n_walkers)
    alive = np.ones(n_walkers, dtype=bool)
    stopped = np.zeros(n_walkers, dtype=bool)
    step = np.zeros(n_walkers, dtype=np.int64)
    stop_step = np.full(n_walkers, -1, dtype=np.int64)
    stop_reason = np.full(n_walkers, -1, dtype=np.int8)
    visited = np.zeros((n_walkers, N), dtype=bool)
    visited[wid, pos_idx] = True
    unique_count = np.ones(n_walkers, dtype=np.int64)
    # aggregate per-site visit counter (summed over walkers)
    visits_agg = np.zeros(N, dtype=np.int64)
    np.add.at(visits_agg, pos_idx, 1)
    # optional per-walker counts (memory: n_walkers * N * 4 bytes)
    visits_pw = (
        np.zeros((n_walkers, N), dtype=np.int32)
        if store_per_walker_visits
        else None
    )
    if visits_pw is not None:
        visits_pw[wid, pos_idx] = 1
    # optional full trajectory history
    history = [pos_idx.copy()] if store_trajectory else None
    # walkers that started already at/above coverage are "stopped-covered" at t=0
    initial = unique_count >= stop_thresh
    stop_step[initial] = 0
    stop_reason[initial] = 0
    stopped |= initial
    return dict(
        wid=wid, alive=alive, stopped=stopped, step=step,
        stop_step=stop_step, stop_reason=stop_reason,
        visited=visited, unique_count=unique_count,
        visits_agg=visits_agg, visits_pw=visits_pw, history=history,
    )


def _finalize_result(
    state: dict,
    pos_idx: np.ndarray,
    neg: np.ndarray,
    frust: np.ndarray,
    Ne: int,
):
    """Pack the walker-state dict into the observable dict returned to callers."""
    return dict(
        unique_visits=state['unique_count'].astype(np.int64),
        unique_frac=state['unique_count'] / neg.shape[0],
        stop_step=state['stop_step'],
        stop_reason=state['stop_reason'],
        final_position=pos_idx.copy(),
        visits_agg=state['visits_agg'],
        visits_per_walker=state['visits_pw'],
        history=state['history'],
        neg_frac=float(neg.sum()) / (2 * Ne),
        frust_frac=float(frust.sum()) / (2 * Ne),
    )


def _absorb_kernel(
    nbrs, neg, frust, *,
    pos_idx, n_walkers, stop_thresh,
    x_node_behavior, rng,
    store_trajectory, store_per_walker_visits,
):
    """Absorbing walker — dies on crossing a killing edge; reflects otherwise."""
    N, z = nbrs.shape
    kill_mask, reflect_mask = _kill_masks(neg, frust, x_node_behavior)
    st = _init_walker_state(
        n_walkers, pos_idx, N, stop_thresh,
        store_trajectory, store_per_walker_visits,
    )
    wid = st['wid']
    # Trap safeguard (reflect mode only): a walker whose position has
    # every outgoing edge reflective AND no killing edge can never move
    # and cannot die — it would infinite-loop. Mark as stopped with
    # reason 2 ('trapped') at step 0.
    if x_node_behavior == 'reflect':
        all_refl = reflect_mask[pos_idx].all(axis=1)
        no_kill = ~kill_mask[pos_idx].any(axis=1)
        trapped = all_refl & no_kill & ~st['stopped']
        if trapped.any():
            st['stop_step'][trapped] = 0
            st['stop_reason'][trapped] = 2
            st['stopped'] |= trapped
    while not st['stopped'].all():
        active = ~st['stopped']
        st['step'] += active
        move = active & st['alive']
        k = rng.integers(0, z, size=n_walkers)
        is_kill = kill_mask[pos_idx, k]
        is_refl = reflect_mask[pos_idx, k]
        # advance only those that neither die-this-step nor reflect-this-step
        advance = move & ~is_refl
        new_pos = np.where(advance, nbrs[pos_idx, k], pos_idx)
        # death is resolved AFTER we know who advanced: absorb-on-crossing
        st['alive'] &= ~(move & is_kill)
        # visit bookkeeping uses the post-move position
        pre = st['visited'][wid, new_pos]
        st['visited'][wid, new_pos] = True
        newly_visited = active & ~pre
        st['unique_count'] += newly_visited
        np.add.at(st['visits_agg'], new_pos[active], 1)
        if st['visits_pw'] is not None:
            st['visits_pw'][wid[active], new_pos[active]] += 1
        pos_idx = new_pos
        if st['history'] is not None:
            st['history'].append(pos_idx.copy())
        _mark_stops(
            st,
            covered=active & (st['unique_count'] >= stop_thresh),
            died=active & ~st['alive'],
        )
    return pos_idx, st


def _kill_kernel(
    nbrs, neg, frust, *,
    pos_idx, n_walkers, stop_thresh,
    x_node_behavior, rng,
    store_trajectory, store_per_walker_visits,
):
    """Killing walker — at each step, die with prob ``1/n_neg`` at the current node."""
    N, z = nbrs.shape
    # kill-count-per-node: under 'reflect' count only frustrated edges,
    # under 'absorb' count all negative edges (notebook-default).
    if x_node_behavior == 'reflect':
        n_kill = frust.sum(axis=1)
    elif x_node_behavior == 'absorb':
        n_kill = neg.sum(axis=1)
    else:
        raise ValueError(f"x_node_behavior must be 'reflect'|'absorb', got {x_node_behavior!r}.")
    kp = np.where(n_kill > 0, 1.0 / np.maximum(n_kill, 1), 0.0)
    st = _init_walker_state(
        n_walkers, pos_idx, N, stop_thresh,
        store_trajectory, store_per_walker_visits,
    )
    wid = st['wid']
    while not st['stopped'].all():
        active = ~st['stopped']
        st['step'] += active
        die = active & st['alive'] & (rng.random(n_walkers) < kp[pos_idx])
        st['alive'] &= ~die
        move = active & st['alive']
        k = rng.integers(0, z, size=n_walkers)
        new_pos = np.where(move, nbrs[pos_idx, k], pos_idx)
        pre = st['visited'][wid, new_pos]
        st['visited'][wid, new_pos] = True
        st['unique_count'] += active & ~pre
        np.add.at(st['visits_agg'], new_pos[active], 1)
        if st['visits_pw'] is not None:
            st['visits_pw'][wid[active], new_pos[active]] += 1
        pos_idx = new_pos
        if st['history'] is not None:
            st['history'].append(pos_idx.copy())
        _mark_stops(
            st,
            covered=active & (st['unique_count'] >= stop_thresh),
            died=active & ~st['alive'],
        )
    return pos_idx, st


def _sticky_kernel(
    nbrs, neg, frust, *,
    pos_idx, n_walkers, stop_thresh,
    x_node_behavior, max_n_cross, rng,
    store_trajectory, store_per_walker_visits,
):
    """Sticky walker — ``P(move) = 2^(-n_cross)``; dies at ``n_cross >= max_n_cross``.

    Real-time is advanced via ``Geometric(2^(-n_cross))`` per iteration
    (event-driven, skipping idle no-move steps).
    """
    N, z = nbrs.shape
    # "crossing" count: under 'reflect' only frustrated counts toward stickiness,
    # under 'absorb' any negative counts.
    if x_node_behavior == 'reflect':
        cross_mask = frust
    elif x_node_behavior == 'absorb':
        cross_mask = neg
    else:
        raise ValueError(f"x_node_behavior must be 'reflect'|'absorb', got {x_node_behavior!r}.")
    st = _init_walker_state(
        n_walkers, pos_idx, N, stop_thresh,
        store_trajectory, store_per_walker_visits,
    )
    wid = st['wid']
    n_cross = np.zeros(n_walkers, dtype=np.int64)
    while not st['stopped'].all():
        active = ~st['stopped']
        p_move = np.ldexp(1.0, -n_cross)
        W = rng.geometric(p_move).astype(np.int64)
        st['step'] += active * W
        k = rng.integers(0, z, size=n_walkers)
        is_cross = cross_mask[pos_idx, k]
        new_pos = np.where(active, nbrs[pos_idx, k], pos_idx)
        n_cross += active & is_cross
        st['alive'] = n_cross < max_n_cross
        pre = st['visited'][wid, new_pos]
        st['visited'][wid, new_pos] = True
        st['unique_count'] += active & ~pre
        np.add.at(st['visits_agg'], new_pos[active], 1)
        if st['visits_pw'] is not None:
            st['visits_pw'][wid[active], new_pos[active]] += 1
        pos_idx = new_pos
        if st['history'] is not None:
            st['history'].append(pos_idx.copy())
        _mark_stops(
            st,
            covered=active & (st['unique_count'] >= stop_thresh),
            died=active & ~st['alive'],
        )
    return pos_idx, st


def _mark_stops(state: dict, covered: np.ndarray, died: np.ndarray) -> None:
    """Record stop time/reason for walkers that just stopped."""
    stopped = state['stopped']
    new_c = covered & ~stopped
    new_d = died & ~stopped
    state['stop_step'][new_c] = state['step'][new_c]
    state['stop_reason'][new_c] = 0
    state['stop_step'][new_d] = state['step'][new_d]
    state['stop_reason'][new_d] = 1
    state['stopped'][new_c | new_d] = True


_KERNELS = {
    'absorb': _absorb_kernel,
    'kill': _kill_kernel,
    'sticky': _sticky_kernel,
}


def run_walker(
    sg: Any,
    rule: str,
    *,
    n_walkers: int,
    start_positions: np.ndarray,
    seed: int,
    coverage_stop: float,
    x_node_behavior: str = 'reflect',
    max_n_cross: int = 14,
    store_trajectory: bool = False,
    store_per_walker_visits: bool = False,
):
    """Run an ensemble of single-walker simulations under one rule.

    Parameters
    ----------
    sg : SignedGraphProtocol
    rule : {'absorb', 'kill', 'sticky'}
    n_walkers : int
    start_positions : (n_walkers,) int — initial node index for each walker
    seed : int
        Seed for a fresh ``np.random.default_rng``; no global state is touched.
    coverage_stop : float
        Walker terminates once it has visited ``coverage_stop * N`` distinct sites.
    x_node_behavior : {'reflect', 'absorb'}
        See module docstring.
    max_n_cross : int
        Sticky-rule death threshold.
    store_trajectory : bool
        If True, ``result['history']`` holds an ``(T+1, n_walkers)`` int64 array.
    store_per_walker_visits : bool
        If True, ``result['visits_per_walker']`` holds ``(n_walkers, N)`` int32 counts.

    Returns
    -------
    dict with keys
        ``unique_visits``, ``unique_frac``, ``stop_step``, ``stop_reason``,
        ``final_position``, ``visits_agg``, ``visits_per_walker``, ``history``,
        ``neg_frac``, ``frust_frac``.
    """
    if rule not in _KERNELS:
        raise ValueError(f"rule must be one of {tuple(_KERNELS)}, got {rule!r}.")
    pos_idx = np.asarray(start_positions, dtype=np.int64).copy()
    if pos_idx.shape != (n_walkers,):
        raise ValueError(
            f"start_positions shape {pos_idx.shape} != ({n_walkers},)."
        )
    nbrs, neg, frust = signed_lattice_tables(sg)
    N = nbrs.shape[0]
    if (pos_idx < 0).any() or (pos_idx >= N).any():
        raise ValueError("start_positions contains out-of-range indices.")
    stop_thresh = max(1, int(coverage_stop * N))
    rng = np.random.default_rng(seed)
    kwargs = dict(
        pos_idx=pos_idx, n_walkers=n_walkers, stop_thresh=stop_thresh,
        x_node_behavior=x_node_behavior, rng=rng,
        store_trajectory=store_trajectory,
        store_per_walker_visits=store_per_walker_visits,
    )
    if rule == 'sticky':
        kwargs['max_n_cross'] = max_n_cross
    final_pos, state = _KERNELS[rule](nbrs, neg, frust, **kwargs)
    Ne = int(sg.Ne) if hasattr(sg, 'Ne') else int(sg.num_edges)
    return _finalize_result(state, final_pos, neg, frust, Ne)


__all__ = [
    'signed_lattice_tables',
    'run_walker',
]
