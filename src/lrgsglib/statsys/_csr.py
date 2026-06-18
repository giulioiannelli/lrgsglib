"""Shared CSR adjacency builders for native (pybind11 / CuPy) statsys backends.

Produces the flat CSR layout consumed by the C++ ``GraphCSR`` helper:

    neigh_indices : int64[total_degree]   concatenated neighbour lists
    neigh_weights : float64[total_degree] concatenated signed edge weights
    neigh_ptr     : int64[N + 1]          row pointers into the two arrays

The undirected graph is symmetrised (each edge contributes both directions).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from ..graphs.protocols import SignedGraphProtocol


def build_graph_csr_from_arrays(
    src: NDArray,
    dst: NDArray,
    wts: NDArray,
    N: int,
) -> tuple[NDArray, NDArray, NDArray]:
    """Build symmetric CSR arrays from directed edge arrays (vectorised)."""
    all_src = np.concatenate([src, dst])
    all_dst = np.concatenate([dst, src])
    all_wts = np.concatenate([wts, wts])

    order = np.argsort(all_src, kind="stable")
    all_dst = all_dst[order]
    all_wts = all_wts[order]
    all_src = all_src[order]

    nodes = np.arange(N + 1, dtype=np.int64)
    ptr = np.searchsorted(all_src, nodes, side="left")

    return all_dst.astype(np.int64), all_wts.astype(np.float64), ptr


def build_graph_csr(
    sg: "SignedGraphProtocol",
    N: int,
) -> tuple[NDArray, NDArray, NDArray]:
    """Build CSR arrays from a signed graph for the native backends.

    Uses the vectorised ``get_edge_arrays()`` path when available (GT graphs)
    and otherwise falls back to the per-node ``get_neighbors_with_weights()``
    loop (fast for NX dict-of-dicts).
    """
    if hasattr(sg, "get_edge_arrays"):
        src, dst, wts = sg.get_edge_arrays()
        return build_graph_csr_from_arrays(src, dst, wts, N)

    indices_list: list[int] = []
    weights_list: list[float] = []
    ptr: list[int] = [0]
    for node in range(N):
        for nn, w in sg.get_neighbors_with_weights(node):
            indices_list.append(nn)
            weights_list.append(w)
        ptr.append(len(indices_list))
    return (
        np.array(indices_list, dtype=np.int64),
        np.array(weights_list, dtype=np.float64),
        np.array(ptr, dtype=np.int64),
    )
