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

import math
from collections import deque
from typing import Any, Callable, List, Optional, Set, Tuple

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


def hub_central_edge(sg, on_g: str) -> Edge:
    """A graph-general central edge: the highest-degree node and an edge to its
    highest-degree neighbour.

    This is the non-geometric fallback for ``get_central_edge`` so the central
    ``single`` / ``singleXERR`` / ``singleZERR`` patterns work on *any* graph
    (ER, BA, Holme-Kim, …), not only lattices (which override with a
    geometry-aware centre). The hub is a natural "centre" for scale-free and
    random graphs. Deterministic given a constructed graph: ties fall to the
    engine's stable node/neighbour iteration order. Raises if the graph has no
    edges.
    """
    nodes = sg.get_nodes_list(on_g)
    best_node, best_deg = None, -1
    for n in nodes:
        d = len(sg.get_graph_neighbors(n, on_g))
        if d > best_deg:
            best_deg, best_node = d, n
    if best_node is None or best_deg == 0:
        raise ValueError(
            "hub_central_edge: graph has no edges; cannot pick a central edge"
        )
    nbrs = sg.get_graph_neighbors(best_node, on_g)
    nb = max(nbrs, key=lambda j: len(sg.get_graph_neighbors(j, on_g)))
    return (best_node, nb)


def star_edges(sg, node: Any, on_g: str) -> List[Edge]:
    """All edges incident to ``node`` — the XERR (star) pattern.

    Returned as ``(node, neighbor)`` tuples (not reordered) to match the
    historical NX ``get_links_XERR`` orientation.
    """
    return [(node, nb) for nb in sg.get_graph_neighbors(node, on_g)]


# --- canonical minimal-cycle ("elementary cell" / ZERR) machinery -----------
#
# The ZERR cell is the smallest cycle (girth face) through ``node``. A regular
# lattice node sits on several equal-length faces, so "the" smallest cycle is
# only defined up to a tie-break. The previous implementation broke ties by
# neighbour-iteration order, which is *engine-dependent*: NX and graph-tool
# visit neighbours in different orders, so they picked different faces (and a
# square-lattice node could even alias to a neighbour's face). We make the
# choice a deterministic function of the node *labels* only — identical on both
# engines — by enumerating every minimal face and selecting canonically.

#: hard caps so the per-seed cost stays ``O(cell)`` even on dense hubs (many
#: equal-length faces / many shortest paths); the deterministic sorted order
#: makes any truncation reproducible across engines.
_MAX_SHORTEST_PATHS = 64
_MAX_FACES = 512


def _path_edges(path: List[Any]) -> List[Edge]:
    """Ordered ``(u, v)`` edges of a vertex ``path`` ``[v0, v1, ...]``."""
    return [_order(path[i], path[i + 1]) for i in range(len(path) - 1)]


def _explore_for_cycles(sg, node: Any, on_g: str):
    """BFS from ``node`` returning ``(dist, preds, L)``.

    ``dist`` maps each reached vertex to its hop distance; ``preds[v]`` is the
    full set of shortest-path predecessors of ``v`` (the BFS shortest-path DAG,
    so *all* shortest paths are recoverable, not just one tree); ``L`` is the
    length of the smallest cycle through ``node`` (or ``None`` if acyclic). The
    BFS early-stops once ``2 * depth > L`` (deeper vertices cannot lie on a
    cycle of length ``<= L``), so the cost is ``O(cell)`` not ``O(graph)``.
    Each vertex also carries an ``arm`` tag — the immediate neighbour of
    ``node`` it descends from — used only to recognise that two predecessors
    reach ``node`` by genuinely distinct routes (a real cycle, not a backtrack).
    """
    dist = {node: 0}
    preds: dict = {node: []}
    arm: dict = {node: None}
    frontier = [node]
    depth = 0
    L = None
    while frontier:
        if L is not None and 2 * depth > L:
            break
        nxt: List[Any] = []
        for u in frontier:
            du = dist[u]
            au = arm[u]
            for w in sorted(sg.get_graph_neighbors(u, on_g)):
                if w == node:
                    continue
                if w not in dist:
                    dist[w] = du + 1
                    preds[w] = [u]
                    # an immediate neighbour of ``node`` is its own arm.
                    arm[w] = w if u == node else au
                    nxt.append(w)
                elif dist[w] == du + 1:
                    # another shortest-path predecessor; distinct arms close an
                    # even cycle of length 2*(du+1) through ``node``.
                    if arm[w] is not None and au is not None and arm[w] != au:
                        cl = 2 * (du + 1)
                        if L is None or cl < L:
                            L = cl
                    preds[w].append(u)
                elif dist[w] == du and au is not None and arm[w] is not None \
                        and arm[w] != au:
                    # same-layer chord between distinct arms → odd cycle 2*du+1.
                    cl = 2 * du + 1
                    if L is None or cl < L:
                        L = cl
        depth += 1
        frontier = nxt
    return dist, preds, L


def _shortest_paths(node: Any, target: Any, preds: dict) -> List[List[Any]]:
    """All shortest paths ``node -> target`` as vertex lists ``[node, ..., target]``,
    rebuilt from the predecessor DAG (bounded by ``_MAX_SHORTEST_PATHS``)."""
    if target == node:
        return [[node]]
    out: List[List[Any]] = []
    stack: List[Tuple[Any, List[Any]]] = [(target, [target])]
    while stack and len(out) < _MAX_SHORTEST_PATHS:
        v, suffix = stack.pop()
        for p in sorted(preds.get(v, ()), reverse=True):
            ext = suffix + [p]
            if p == node:
                out.append(list(reversed(ext)))
                if len(out) >= _MAX_SHORTEST_PATHS:
                    break
            else:
                stack.append((p, ext))
    return out


def minimal_cycles_through(sg, node: Any, on_g: str) -> List[List[Edge]]:
    """Every minimal cycle (girth face) through ``node``, each as a sorted edge
    list. Deterministic and engine-independent (depends only on node labels and
    adjacency). Returns ``[]`` when ``node`` lies on no cycle (tree / pendant).

    This is the shared primitive behind both the geometry-free canonical ZERR
    cell (:func:`elementary_cell_edges`, the lexicographically-smallest face)
    and the lattice oriented-plaquette picker (which selects among these faces
    by coordinate orientation).
    """
    dist, preds, L = _explore_for_cycles(sg, node, on_g)
    if L is None:
        return []
    faces: dict = {}  # frozenset(edges) -> sorted edge list

    def _add(edges: List[Edge]) -> None:
        key = frozenset(edges)
        if len(key) == L and key not in faces:
            faces[key] = sorted(key)

    if L % 2 == 0:
        r = L // 2
        for w in sorted(v for v, d in dist.items() if d == r):
            paths = _shortest_paths(node, w, preds)
            for i in range(len(paths)):
                si = set(paths[i])
                for j in range(i + 1, len(paths)):
                    # two shortest node->w paths sharing only the endpoints
                    # bound a 2r-cycle.
                    if si & set(paths[j]) == {node, w}:
                        _add(_path_edges(paths[i]) + _path_edges(paths[j]))
                        if len(faces) >= _MAX_FACES:
                            return sorted(faces.values())
    else:
        r = (L - 1) // 2
        layer = sorted(v for v, d in dist.items() if d == r)
        layerset = set(layer)
        seen: Set[Edge] = set()
        for u in layer:
            for w in sorted(sg.get_graph_neighbors(u, on_g)):
                if w not in layerset or w == u:
                    continue
                chord = _order(u, w)
                if chord in seen:
                    continue
                seen.add(chord)
                pu = _shortest_paths(node, u, preds)
                pw = _shortest_paths(node, w, preds)
                for pa in pu:
                    sa = set(pa)
                    for pb in pw:
                        # disjoint arms (share only ``node``) → 2r+1 cycle.
                        if sa & set(pb) == {node}:
                            _add(_path_edges(pa) + _path_edges(pb) + [chord])
                            if len(faces) >= _MAX_FACES:
                                return sorted(faces.values())
    return sorted(faces.values())


def elementary_cell_edges(sg, node: Any, on_g: str) -> List[Edge]:
    """Canonical edges of the smallest cycle through ``node`` — the generalized
    ZERR cell (4-cycle on a square lattice, 3-cycle on triangular/kagome,
    6-cycle on honeycomb, …).

    Among all minimal faces through ``node`` the **lexicographically-smallest**
    one (compared as a sorted edge-tuple list) is returned, so the choice is a
    deterministic function of the node labels and therefore *identical on the
    NetworkX and graph-tool engines*. Returns ``[]`` if ``node`` lies on no
    cycle. This is the geometry-free baseline; lattices refine it with a
    coordinate-oriented picker (see ``SignedGraph.cell_edges``) that additionally
    makes distinct seed nodes own distinct faces.
    """
    faces = minimal_cycles_through(sg, node, on_g)
    return min(faces) if faces else []


def _min_image(d: float, period: Optional[float]) -> float:
    """Nearest-image displacement: fold ``d`` into ``[-period/2, period/2]``."""
    if not period:
        return d
    return d - period * round(d / period)


def oriented_cell_edges(
    sg,
    node: Any,
    on_g: str,
    posfn: Callable[[Any], Optional[Tuple[float, float]]],
    box: Optional[Tuple[Optional[float], Optional[float]]] = None,
) -> List[Edge]:
    """Coordinate-oriented ZERR cell — the lattice refinement of
    :func:`elementary_cell_edges`.

    Among the minimal faces through ``node`` (the same set the geometry-free
    baseline chooses its lex-min from), pick the one whose centroid lies at the
    smallest counter-clockwise angle from the ``+x`` axis as seen from ``node``.
    This is a *fixed rotation*: on a regular lattice every node owns the same
    "first" plaquette (e.g. its ``+x/+y`` / north-east face on a square lattice),
    so distinct seed nodes own distinct cells in the bulk — which the pure
    topological lex-min cannot guarantee (adjacent nodes can share the same
    lex-min face). ``posfn`` maps a node (in the ``on_g`` representation) to its
    ``(x, y)`` coordinate; where it returns ``None`` (no geometry / boundary)
    the function degrades to the canonical lex-min face, so it is always safe to
    call. ``box`` gives the periods ``(Lx, Ly)`` for periodic lattices so the
    angle is measured on the *nearest image* of each face vertex — without it a
    plaquette straddling the periodic seam would point the wrong way and alias
    onto a neighbour's cell. Ties in angle fall back to the lex-min face for
    full determinism.
    """
    faces = minimal_cycles_through(sg, node, on_g)
    if not faces:
        return []
    if len(faces) == 1:
        return faces[0]
    p0 = posfn(node)
    if p0 is None:
        return min(faces)
    bx, by = box if box is not None else (None, None)

    def _angle_key(face: List[Edge]):
        verts = {v for e in face for v in e}
        cx = cy = 0.0
        for v in verts:
            p = posfn(v)
            if p is None:
                return (math.inf, face)
            cx += _min_image(p[0] - p0[0], bx)
            cy += _min_image(p[1] - p0[1], by)
        cx /= len(verts)
        cy /= len(verts)
        ang = math.atan2(cy, cx) % (2.0 * math.pi)
        return (ang, face)

    return min(faces, key=_angle_key)


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
