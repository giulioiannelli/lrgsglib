"""Shared CSR adjacency builders for native (pybind11 / CuPy) statsys backends.

Produces the flat CSR layout consumed by the C++ ``GraphCSR`` helper:

    neigh_indices : int64[total_degree]   concatenated neighbour lists
    neigh_weights : float64[total_degree] concatenated signed edge weights
    neigh_ptr     : int64[N + 1]          row pointers into the two arrays

The undirected graph is symmetrised (each edge contributes both directions).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from ..graphs.protocols import SignedGraphProtocol

#: How a spin model's symmetric coupling matrix J derives from the signed
#: edge weights (Ising plan §3b, shared by every spin substrate):
#: 'raw' J_ij = w_ij; 'sym' J_ij = w_ij/sqrt(deg_i deg_j) (normalized
#: adjacency); 'avg' J_ij = w_ij (1/deg_i + 1/deg_j) (per-site-averaged,
#: self-symmetrized).
COUPLING_NORMS: tuple[str, ...] = ("raw", "sym", "avg")


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


class CouplingCSR(NamedTuple):
    """The ONE symmetric coupling structure of a spin model (plan §3b):
    the CSR adjacency with the coupling normalization already applied — the
    same J feeds ΔE and the energy observable by construction — plus the
    row-index expansion and each undirected edge listed once (u < v) for
    whole-bond sweeps (Swendsen–Wang)."""

    idx: NDArray  #: int64 concatenated neighbour lists
    J: NDArray  #: float64 normalized couplings, aligned with ``idx``
    ptr: NDArray  #: int64 row pointers (N + 1)
    rows: NDArray  #: int64 row index of each CSR entry
    edge_u: NDArray  #: int64 undirected-edge endpoints (u < v)
    edge_v: NDArray
    edge_J: NDArray  #: float64 coupling of each undirected edge


def build_coupling_csr(
    sg: "SignedGraphProtocol",
    N: int,
    norm: str = "raw",
) -> CouplingCSR:
    """Build a spin model's symmetric coupling CSR (see :data:`COUPLING_NORMS`).

    Shared home for the coupling construction every spin substrate
    (``IsingBase``, ``PottsBase``, ...) performs once at init.
    """
    if norm not in COUPLING_NORMS:
        raise ValueError(
            f"coupling norm must be one of {COUPLING_NORMS}; got {norm!r}."
        )
    idx, wts, ptr = build_graph_csr(sg, N)
    deg = np.diff(ptr)
    rows = np.repeat(np.arange(N, dtype=np.int64), deg)
    if norm == "raw":
        J = wts
    elif norm == "sym":
        J = wts / np.sqrt(deg[rows] * deg[idx])
    else:  # 'avg' — per-site-averaged, self-symmetrized
        J = wts * (1.0 / deg[rows] + 1.0 / deg[idx])
    J = np.asarray(J, dtype=np.float64)
    upper = rows < idx
    return CouplingCSR(
        idx=idx,
        J=J,
        ptr=ptr,
        rows=rows,
        edge_u=rows[upper],
        edge_v=idx[upper],
        edge_J=J[upper],
    )
