"""Edge generation utilities for GraphOfGraphs.

These functions generate edge lists for the composite graph structure,
handling base, fiber, and anchor edges separately.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    import graph_tool as gt


def extract_edges_from_gt_graph(
    graph: "gt.Graph",
    vertex_offset: int = 0,
) -> list[tuple[int, int]]:
    """Extract edge list from a graph-tool graph with optional vertex offset.

    Parameters
    ----------
    graph : gt.Graph
        Source graph to extract edges from
    vertex_offset : int, default 0
        Offset to add to all vertex indices

    Returns
    -------
    list[tuple[int, int]]
        Edge list as (source, target) tuples with offset applied
    """
    return [
        (int(e.source()) + vertex_offset, int(e.target()) + vertex_offset)
        for e in graph.edges()
    ]


def generate_fiber_edges(
    fiber_graph: "gt.Graph",
    base_idx: int,
    n_base: int,
    n_fiber: int,
) -> list[tuple[int, int]]:
    """Generate edges for a fiber attached to a specific base node.

    Parameters
    ----------
    fiber_graph : gt.Graph
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
    return extract_edges_from_gt_graph(fiber_graph, vertex_offset=fiber_start)


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


def collect_edge_signs(
    graph: "gt.Graph",
    edge_list: list[tuple[int, int]],
    vertex_offset: int = 0,
) -> list[int]:
    """Extract edge signs from a graph-tool graph for a given edge list.

    Parameters
    ----------
    graph : gt.Graph
        Source graph with 'sign' edge property
    edge_list : list[tuple[int, int]]
        Edges to get signs for (with offset already applied)
    vertex_offset : int, default 0
        Offset that was applied to edges (to reverse for lookup)

    Returns
    -------
    list[int]
        Sign values (+1 or -1) for each edge
    """
    if "sign" not in graph.edge_properties:
        return [1] * len(edge_list)

    sign_prop = graph.edge_properties["sign"]
    signs = []

    for u, v in edge_list:
        # Reverse offset for lookup
        u_local = u - vertex_offset
        v_local = v - vertex_offset
        e = graph.edge(u_local, v_local)
        if e is not None:
            signs.append(sign_prop[e])
        else:
            signs.append(1)

    return signs
