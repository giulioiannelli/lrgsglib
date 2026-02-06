"""Basic graph manipulation utilities."""

from __future__ import annotations

import warnings
from typing import List, Tuple

import numpy as np
from networkx import (
    Graph,
    connected_components,
    is_connected,
    selfloop_edges,
    set_edge_attributes,
)

__all__ = [
    "fast_set_weights_from_matrix",
    "remove_edges",
    "rewire_edges_optimized",
    "get_giant_component",
    "get_giant_component_leftoff",
]


def fast_set_weights_from_matrix(G: Graph, weight_matrix: np.ndarray) -> None:
    """Update edge weights of a fully connected graph from a matrix."""
    num_nodes = G.number_of_nodes()
    if weight_matrix.shape != (num_nodes, num_nodes):
        raise ValueError(
            "The dimensions of weight_matrix must match the number of nodes in G."
        )

    edge_weights = {
        (i, j): weight_matrix[i, j]
        for i in range(num_nodes)
        for j in range(i + 1, num_nodes)
    }
    set_edge_attributes(G, edge_weights, "weight")


def remove_edges(G: Graph, pdil: float) -> Graph:
    """Remove a fraction of edges while keeping the graph connected."""
    edges = list(G.edges())
    num_edges_to_remove = int(len(edges) * pdil)
    rndE = np.random.choice(len(edges), size=num_edges_to_remove, replace=False)
    edges_to_remove = [edges[i] for i in rndE]
    G.remove_edges_from(edges_to_remove)

    if not is_connected(G):
        warnings.warn(
            "The resulting graph is no longer connected. Returning the largest connected component."
        )
        giant_component = max(connected_components(G), key=len)
        G = G.subgraph(giant_component).copy()
        return G

    return G


def rewire_edges_optimized(G: Graph, prew: float) -> Graph:
    """Rewire edges with probability ``prew`` without introducing duplicates."""
    nodes = list(G.nodes())
    edges = list(G.edges())
    num_edges = len(edges)
    rewire_flags = np.random.rand(num_edges) < prew
    existing_edges = set(G.edges())

    for i, edge in enumerate(edges):
        if rewire_flags[i]:
            u, v = edge
            G.remove_edge(u, v)
            existing_edges.remove((u, v))

            new_neighbor = nodes[np.random.randint(0, len(nodes))]
            attempts = 0
            while (
                new_neighbor == u
                or new_neighbor == v
                or (u, new_neighbor) in existing_edges
                or (new_neighbor, u) in existing_edges
            ):
                new_neighbor = nodes[np.random.randint(0, len(nodes))]
                attempts += 1
                if attempts > len(nodes):
                    break

            if new_neighbor != u and not G.has_edge(u, new_neighbor):
                G.add_edge(u, new_neighbor)
                existing_edges.add((u, new_neighbor))
    return G


def get_giant_component(graph: Graph, remove_selfloops: bool = True) -> Graph:
    """Return the largest connected component of ``graph``."""
    G_clean = graph.copy()
    if remove_selfloops:
        G_clean.remove_edges_from(selfloop_edges(G_clean))
    components = sorted(connected_components(G_clean), key=len, reverse=True)
    return G_clean.subgraph(components[0]).copy()


def get_giant_component_leftoff(
    graph: Graph, remove_selfloops: bool = True
) -> Tuple[Graph, List]:
    """Return the giant component and nodes left outside of it."""
    G_clean = graph.copy()
    if remove_selfloops:
        G_clean.remove_edges_from(selfloop_edges(G_clean))

    components = sorted(connected_components(G_clean), key=len, reverse=True)
    if components:
        giant_nodes = components[0]
    else:
        giant_nodes = set()

    giant_component = G_clean.subgraph(giant_nodes).copy()
    left_off_nodes = list(set(G_clean.nodes()) - giant_nodes)

    return giant_component, left_off_nodes
