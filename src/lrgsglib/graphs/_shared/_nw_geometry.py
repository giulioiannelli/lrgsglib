"""Engine-neutral geometry helpers for negative-link (``nwDict``) patterns.

These functions express the *geometry* of the defect patterns (star / XERR,
elementary cell / ZERR, radius-R ball) in terms of a single cross-engine seam:
``sg.get_graph_neighbors(node, on_g)``. They therefore run unchanged on both the
NetworkX and graph-tool engines, where ``sg`` is a ``SignedGraphNX`` /
``SignedGraphGT`` and nodes are whatever that engine uses for the ``on_g``
representation (integers for GT and for the NX ``'G'`` representation).

All edge-returning helpers yield **sorted ``(u, v)`` tuples** so results are
order-independent and directly comparable across engines.

Design notes
------------
* ``elementary_cell_edges`` is the generalized ZERR primitive: it returns the
  edges of the *smallest cycle through* ``node`` (its girth face). On a square
  lattice this is a 4-cycle, on a triangular/kagome lattice a 3-cycle, on a
  honeycomb a 6-cycle, etc. — so a single helper covers every geometry with a
  uniform cell size and degrades gracefully on mixed-cell lattices (it always
  picks the smallest face incident to the node). This mirrors the intent of
  ``lrgsglib.graphs.nx.funcs.neighbors.get_smallest_cycle_graph_node`` but is
  engine-neutral.
* ``neighbors_at_distance`` matches ``get_neighbors_at_distance`` exactly:
  nodes at *exactly* ``R`` hops (not the closed ball), which is what the NX
  ``get_links_rball`` consumes.
"""
from __future__ import annotations

from collections import deque
from typing import Any, List, Set, Tuple

# A representation key is engine-defined; default resolution is left to the
# caller (containers pass an explicit ``on_g``).

Edge = Tuple[Any, Any]


def _order(u: Any, v: Any) -> Edge:
    """Return ``(u, v)`` ordered so the tuple is canonical for set membership."""
    return (u, v) if u <= v else (v, u)


def common_neighbors(sg, a: Any, b: Any, on_g: str) -> Set[Any]:
    """Nodes adjacent to both ``a`` and ``b`` (engine-neutral).

    Replaces ``networkx.common_neighbors``.
    """
    na = set(sg.get_graph_neighbors(a, on_g))
    nb = set(sg.get_graph_neighbors(b, on_g))
    return (na & nb) - {a, b}


def induced_edges(sg, nodes, on_g: str) -> List[Edge]:
    """Edges of the subgraph induced by ``nodes`` (sorted ``(u, v)`` tuples).

    Replaces ``graph.subgraph(nodes).edges()``.
    """
    nodeset = set(nodes)
    seen: Set[Edge] = set()
    for u in nodeset:
        for w in sg.get_graph_neighbors(u, on_g):
            if w in nodeset and u != w:
                seen.add(_order(u, w))
    return sorted(seen)


def neighbors_at_distance(sg, node: Any, R: int, on_g: str) -> List[Any]:
    """Nodes at *exactly* ``R`` hops from ``node`` (BFS).

    Matches ``lrgsglib.graphs.nx.funcs.neighbors.get_neighbors_at_distance``.
    """
    dist = {node: 0}
    q = deque([node])
    while q:
        u = q.popleft()
        for w in sg.get_graph_neighbors(u, on_g):
            if w not in dist:
                dist[w] = dist[u] + 1
                q.append(w)
    return [n for n, d in dist.items() if d == R]


def star_edges(sg, node: Any, on_g: str) -> List[Edge]:
    """All edges incident to ``node`` — the XERR (star) pattern.

    Returned as ``(node, neighbor)`` tuples (not reordered) to match the
    historical NX ``get_links_XERR`` orientation.
    """
    return [(node, nb) for nb in sg.get_graph_neighbors(node, on_g)]


def elementary_cell_edges(sg, node: Any, on_g: str) -> List[Edge]:
    """Edges of the smallest cycle through ``node`` — the generalized ZERR cell.

    Layered BFS outward from ``node`` recording, for every reached vertex, which
    *arm* (which immediate neighbour of ``node``) discovered it. A non-tree edge
    joining two different arms closes a cycle through ``node`` of length
    ``dist[u] + dist[w] + 1``; the globally shortest such cycle is the
    elementary face. The BFS stops as soon as no shorter cycle can still appear
    (``2 * depth >= best_length``), so each call costs ``O(cell)`` rather than
    ``O(graph)`` — essential because the ``rand*`` patterns call this once per
    seed node. Returns the cycle's edges as sorted ``(u, v)`` tuples, or ``[]``
    if ``node`` lies on no cycle (e.g. a tree / pendant).
    """
    dist = {node: 0}
    parent: dict = {node: None}
    arm: dict = {node: None}
    best = None  # (length, u, w)
    processed: Set[Edge] = set()
    frontier = [node]
    depth = 0
    while frontier and (best is None or 2 * depth < best[0]):
        nxt: List[Any] = []
        for u in frontier:
            for w in sg.get_graph_neighbors(u, on_g):
                key = _order(u, w)
                if w not in dist:
                    dist[w] = depth + 1
                    parent[w] = u
                    arm[w] = w if u == node else arm[u]
                    nxt.append(w)
                elif (
                    key not in processed
                    and u != node
                    and w != node
                    and arm[u] is not None
                    and arm[w] is not None
                    and arm[u] != arm[w]
                ):
                    length = dist[u] + dist[w] + 1
                    if best is None or length < best[0]:
                        best = (length, u, w)
                processed.add(key)
        depth += 1
        frontier = nxt

    if best is None:
        return []

    _, x, y = best

    def path_to_root(t: Any) -> List[Any]:
        p = []
        while t is not None:
            p.append(t)
            t = parent[t]
        return p  # [t, ..., node]

    arm_x = list(reversed(path_to_root(x)))  # [node, ..., x]
    arm_y = path_to_root(y)                  # [y, ..., node]

    edges: List[Edge] = []
    for i in range(len(arm_x) - 1):
        edges.append(_order(arm_x[i], arm_x[i + 1]))
    edges.append(_order(x, y))
    for i in range(len(arm_y) - 1):
        edges.append(_order(arm_y[i], arm_y[i + 1]))
    # dedup, preserve order
    return list(dict.fromkeys(edges))


def ball_edges(sg, center: Any, R: int, on_g: str) -> Set[Edge]:
    """All edges incident to nodes at exactly distance ``R`` from ``center``.

    Engine-neutral equivalent of the NX ``get_links_rball`` body.
    """
    shell = neighbors_at_distance(sg, center, R, on_g)
    return {
        _order(node, nb)
        for node in shell
        for nb in sg.get_graph_neighbors(node, on_g)
    }
