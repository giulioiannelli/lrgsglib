"""Cluster-size-distribution observable for the voter model.

A *cluster* (domain) is a connected component of the subgraph of **active**
edges. The activation predicate is ``b_ij * s_i * s_j > 0`` where the per-edge
sign ``b_ij`` selects the cluster mode (see ``defaults.py``):

* ``satisfied`` -- ``b_ij = sign(w_ij)``: domains of aligned signed opinion (a
  negative edge "agrees" when the spins are opposite). Reduces to same-opinion
  domains on an unsigned graph.
* ``rawspin``   -- ``b_ij = +1``: domains of equal raw spin value.

The distribution is recomputed from scratch at each recorded sweep (a single
connected-components pass over the active-edge subgraph, O(N + E)). Because the
voter records once per sweep and performs ~N single-spin events per sweep, a
per-record recompute is the same asymptotic cost as maintaining the partition
incrementally -- O(E) per recorded sweep -- but is trivially correct and free of
the fragile, super-linear split bookkeeping an incremental scheme needs under
deletions. The native C/pybind kernel (`_ccore/LRGSG_clusters.{c,h}`) mirrors
this exact computation and is checked against these functions for parity.
"""

from __future__ import annotations

from collections import Counter

import numpy as np


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


def cluster_components(s, idx, b) -> np.ndarray:
    """Component label per node via BFS over active edges.

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
    """Counter ``{cluster_size: number_of_clusters}``."""
    label = cluster_components(s, idx, b)
    sizes = np.bincount(label)
    return Counter(int(x) for x in sizes)
