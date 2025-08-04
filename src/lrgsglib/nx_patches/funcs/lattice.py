"""Lattice graph generation utilities."""

from __future__ import annotations

from typing import Tuple

import numpy as np
from networkx import Graph

from ...utils.basic.iterables import cProd_Iter

__all__ = [
    "LatticeND_graph_FastPatch",
    "LatticeND_graph_with_dilution",
]


def LatticeND_graph_FastPatch(dim: Tuple[int, ...], periodic: bool = False) -> Graph:
    """Generate an N-dimensional cubic lattice graph."""
    if len(dim) == 2 and 2 in dim:
        raise ValueError(
            "A dimension of size 2 in a 2D grid can cause incorrect periodic edge handling."
        )

    G = Graph()
    num_dimensions = len(dim)
    nodes = list(cProd_Iter(dim))
    G.add_nodes_from(nodes)

    e_i = [tuple(1 if i == j else 0 for j in range(num_dimensions)) for i in range(num_dimensions)]

    for pt in nodes:
        for drt in e_i:
            neighbor = tuple((d + p) for d, p in zip(pt, drt))
            if all(0 <= n < dim[i] for i, n in enumerate(neighbor)):
                G.add_edge(pt, neighbor)
            elif periodic:
                neighbor = tuple((n % dim[i]) for i, n in enumerate(neighbor))
                G.add_edge(pt, neighbor)

    return G


def LatticeND_graph_with_dilution(
    dim: Tuple[int, ...], periodic: bool = False, pdil: float = 0.0
) -> Graph:
    """Generate a diluted N-dimensional cubic lattice graph."""
    G = LatticeND_graph_FastPatch(dim, periodic)

    edges = list(G.edges())
    num_edges = len(edges)
    num_edges_to_remove = int(pdil * num_edges)
    edges_to_remove = np.random.choice(num_edges, size=num_edges_to_remove, replace=False)

    for edge_idx in edges_to_remove:
        G.remove_edge(*edges[edge_idx])

    return G
