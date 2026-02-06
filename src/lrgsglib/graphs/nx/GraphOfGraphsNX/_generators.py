"""Edge generation utilities for GraphOfGraphs (NetworkX version).

These functions generate edges for the composite graph structure,
handling base, fiber, and anchor edges separately.
"""

from __future__ import annotations

from typing import Sequence

import networkx as nx


def extract_edges_from_nx_graph(
    graph: nx.Graph,
    vertex_offset: int = 0,
) -> list[tuple[int, int]]:
    """Extract edge list from a NetworkX graph with optional vertex offset.

    Parameters
    ----------
    graph : nx.Graph
        Source graph to extract edges from (assumes integer node labels)
    vertex_offset : int, default 0
        Offset to add to all vertex indices

    Returns
    -------
    list[tuple[int, int]]
        Edge list as (source, target) tuples with offset applied
    """
    edges = []
    for u, v in graph.edges():
        edges.append((int(u) + vertex_offset, int(v) + vertex_offset))
    return edges


def extract_edge_signs_from_nx_graph(
    graph: nx.Graph,
) -> dict[tuple[int, int], int]:
    """Extract edge signs from a NetworkX graph.

    Parameters
    ----------
    graph : nx.Graph
        Source graph with optional 'sign' edge attribute

    Returns
    -------
    dict[tuple[int, int], int]
        Mapping of (u, v) -> sign for all edges
    """
    signs = {}
    for u, v, data in graph.edges(data=True):
        sign = data.get("sign", 1)
        signs[(int(u), int(v))] = sign
    return signs


def generate_fiber_edges(
    fiber_graph: nx.Graph,
    base_idx: int,
    n_base: int,
    n_fiber: int,
) -> list[tuple[int, int]]:
    """Generate edges for a fiber attached to a specific base node.

    Parameters
    ----------
    fiber_graph : nx.Graph
        Prototype fiber graph to copy edges from
    base_idx : int
        Index of the base node this fiber attaches to
    n_base : int
        Total number of base nodes
    n_fiber : int
        Number of nodes in each fiber

    Returns
    -------
    list[tuple[int, int]]
        Edge list for this fiber with global vertex indices
    """
    # Fiber vertex offset: N_base + base_idx * N_fiber
    fiber_start = n_base + base_idx * n_fiber
    return extract_edges_from_nx_graph(fiber_graph, vertex_offset=fiber_start)


def generate_anchor_edges(
    n_base: int,
    n_fiber: int,
    anchor_indices: Sequence[int],
) -> list[tuple[int, int]]:
    """Generate edges connecting fiber anchors to base nodes.

    Parameters
    ----------
    n_base : int
        Number of base nodes
    n_fiber : int
        Number of nodes in each fiber
    anchor_indices : Sequence[int]
        Anchor index within each fiber (one per base node)

    Returns
    -------
    list[tuple[int, int]]
        Anchor edges connecting base nodes to fiber anchors
    """
    edges = []
    for base_idx in range(n_base):
        fiber_start = n_base + base_idx * n_fiber
        anchor_vertex = fiber_start + anchor_indices[base_idx]
        edges.append((base_idx, anchor_vertex))
    return edges


def build_composite_nx_graph(
    n_total: int,
    base_edges: list[tuple[int, int]],
    fiber_edges_list: list[list[tuple[int, int]]],
    anchor_edges: list[tuple[int, int]],
    base_signs: dict[tuple[int, int], int],
    fiber_signs_list: list[dict[tuple[int, int], int]],
    n_base: int,
    n_fiber: int,
) -> nx.Graph:
    """Build the composite NetworkX graph with all edges.

    Parameters
    ----------
    n_total : int
        Total number of nodes
    base_edges : list[tuple[int, int]]
        Edges from the base graph
    fiber_edges_list : list[list[tuple[int, int]]]
        Edges from each fiber (one list per base node)
    anchor_edges : list[tuple[int, int]]
        Anchor edges connecting fibers to base
    base_signs : dict
        Signs for base edges
    fiber_signs_list : list[dict]
        Signs for each fiber's edges
    n_base : int
        Number of base nodes
    n_fiber : int
        Number of fiber nodes

    Returns
    -------
    nx.Graph
        Composite graph with 'sign' edge attribute
    """
    G = nx.Graph()
    G.add_nodes_from(range(n_total))

    # Add base edges
    for u, v in base_edges:
        sign = base_signs.get((u, v), base_signs.get((v, u), 1))
        G.add_edge(u, v, sign=sign)

    # Add fiber edges
    for base_idx, fiber_edges in enumerate(fiber_edges_list):
        fiber_start = n_base + base_idx * n_fiber
        fiber_signs = fiber_signs_list[base_idx]
        for u, v in fiber_edges:
            # Map back to local indices for sign lookup
            u_local = u - fiber_start
            v_local = v - fiber_start
            sign = fiber_signs.get((u_local, v_local), fiber_signs.get((v_local, u_local), 1))
            G.add_edge(u, v, sign=sign)

    # Add anchor edges (always positive)
    for u, v in anchor_edges:
        G.add_edge(u, v, sign=1)

    return G
