"""
Native coordinate generators for :class:`Lattice2DGT`.

Pure-Python, **NetworkX-free** reimplementations of the exotic 2D lattice
geometries (rhomb-octagonal, kagome, tri-hexagonal) so the graph-tool engine is
self-contained.  Each generator returns ``(n_nodes, edges, pos)`` where

* ``n_nodes`` : int, number of integer-labelled vertices ``0 .. n_nodes-1``;
* ``edges``   : list of ``(int, int)`` endpoint pairs;
* ``pos``     : ``(n_nodes, 2)`` float array of plotting coordinates.

The produced topology matches the corresponding NetworkX generator in
``lrgsglib.graphs.nx.Lattice2DNX.generators_2d`` *exactly* (verified by adjacency
spectrum), so the ``nx`` and ``gt`` engines are interchangeable through the
unified :class:`lrgsglib.graphs.Lattice2D` factory.

The kagome and tri-hexagonal lattices are defined in the NX code as the
line-graph / vertex-truncation of NetworkX's *built-in*
``networkx.hexagonal_lattice_graph``; :func:`_builtin_hex` reproduces that exact
node/edge set combinatorially so the derived lattices stay bit-for-bit
compatible without importing NetworkX here.
"""

from __future__ import annotations

from math import sqrt
from typing import Dict, List, Tuple

import numpy as np

# --- geometric layout constants (drawing only; mirror generators_2d.py) -------
_RHOMB_SIZE = 0.4  # half-diagonal of an oct_sqr rhomb
_RHOMB_SPACING = 4.0 * _RHOMB_SIZE
_TRIHEX_SHRINK = 0.25  # truncation offset toward each hex neighbour
_HEX_H = sqrt(3.0) / 2.0  # row height of the hexagonal embedding

Coord = Tuple[int, int]
GenResult = Tuple[int, List[Tuple[int, int]], np.ndarray]


def _finalize(nodes, edges_coord, pos_fn) -> GenResult:
    """Relabel coordinate nodes to ``0..N-1`` (sorted), dedup edges, build pos.

    Self-loops and duplicate edges (which periodic wrapping can produce) are
    dropped, mirroring ``networkx.Graph`` semantics.
    """
    order = sorted(nodes)
    idx = {c: i for i, c in enumerate(order)}
    eset = set()
    for a, b in edges_coord:
        ia, ib = idx[a], idx[b]
        if ia == ib:
            continue
        eset.add((ia, ib) if ia < ib else (ib, ia))
    pos = np.array([pos_fn(c) for c in order], dtype=float)
    return len(order), sorted(eset), pos


# --- squared -----------------------------------------------------------------
def squared(m: int, n: int, periodic: bool = False) -> GenResult:
    """Square lattice (4-neighbour grid); bulk coordination ``z = 4``.

    Reproduces ``...nx...generators_2d.squared_lattice_graph_FastPatch``.
    """
    nodes = {(i, j) for i in range(m) for j in range(n)}
    ec: List[Tuple[Coord, Coord]] = []
    for i in range(m - 1):
        for j in range(n):
            ec.append(((i, j), (i + 1, j)))
    for i in range(m):
        for j in range(n - 1):
            ec.append(((i, j), (i, j + 1)))
    if periodic and m > 2:
        ec += [((0, j), (m - 1, j)) for j in range(n)]
    if periodic and n > 2:
        ec += [((i, 0), (i, n - 1)) for i in range(m)]
    return _finalize(nodes, ec, lambda c: (float(c[0]), float(c[1])))


# --- triangular --------------------------------------------------------------
def triangular(m: int, n: int, periodic: bool = False) -> GenResult:
    """Triangular lattice; bulk coordination ``z = 6``.

    Reproduces ``...nx...generators_2d.triangular_lattice_graph_FastPatch``;
    periodic boundaries require ``m >= 3`` and ``n >= 5`` (as in NX).
    """
    if periodic and (m < 3 or n < 5):
        raise ValueError(
            f"periodic triangular needs side1 >= 3 and side2 >= 5, "
            f"got side1={m}, side2={n}"
        )
    nodes: set = set()
    ec: List[Tuple[Coord, Coord]] = []

    def add(a: Coord, b: Coord) -> None:
        ec.append((a, b))
        nodes.add(a)
        nodes.add(b)

    if periodic:
        for i in range(n + 1):
            for j in range(m + 1):
                nodes.add((i % n, j % m))
                if i < n or j % 2 == 0:
                    add((i % n, j % m), ((i + 1) % n, j % m))
                if j < m:
                    add((i % n, j % m), (i % n, (j + 1) % m))
                    if j % 2:
                        add((i % n, j % m), ((i + 1) % n, (j + 1) % m))
                    else:
                        add((i % n, j % m), ((i - 1) % n, (j + 1) % m))
    else:
        nodes = {(i, j) for i in range(n) for j in range(m)}
        for j in range(m):
            for i in range(n - 1):
                add((i, j), (i + 1, j))
        for j in range(m - 1):
            for i in range(n):
                add((i, j), (i, j + 1))
        for j in range(1, m - 1, 2):
            for i in range(n - 1):
                add((i, j), (i + 1, j + 1))
        for j in range(0, m - 1, 2):
            for i in range(n - 1):
                add((i + 1, j), (i, j + 1))

    return _finalize(
        nodes, ec, lambda c: (0.5 * (c[1] % 2) + c[0], _HEX_H * c[1])
    )


# --- hexagonal ---------------------------------------------------------------
def hexagonal(n: int, m: int, periodic: bool = False) -> GenResult:
    """Hexagonal (honeycomb) lattice; bulk coordination ``z = 3``.

    Reproduces ``...nx...generators_2d.hexagonal_lattice_graph_FastPatch``.
    Note the NX argument order ``(n, m) = (side1, side2)``.
    """
    nodes: set = set()
    ec: List[Tuple[Coord, Coord]] = []

    def add(a: Coord, b: Coord) -> None:
        ec.append((a, b))
        nodes.add(a)
        nodes.add(b)

    for i in range(n):
        for j in range(m - 1):
            add((i, j), (i, j + 1))
    for i in range(n - 1):
        for j in range(m):
            if i % 2 == j % 2:
                add((i, j), (i + 1, j))
    if periodic:
        for j in range(m):
            if j % 2:
                add((0, j), (n - 1, j))
        for i in range(n):
            if (i + 1) % 2:
                add((i, 0), (i, m - 1))

    def pos_fn(c: Coord) -> Tuple[float, float]:
        i, j = c
        return (0.5 + i + i // 2 + (j % 2) * ((i % 2) - 0.5), _HEX_H * j)

    return _finalize(nodes, ec, pos_fn)


# --- rhomb-octagonal (octagonal_sqr) -----------------------------------------
def rhomb_octagonal(m: int, n: int, periodic: bool = False) -> GenResult:
    """Rhomb-octagonal lattice: every square-lattice site becomes a 4-node rhomb.

    ``N = 4 * m * n``; bulk coordination ``z = 3``.  Integer node layout per
    rhomb ``base = 4*idx``: ``0`` bottom, ``1`` left, ``2`` right, ``3`` top,
    matching :func:`...nx...generators_2d.rhomb_octagonal_graph_FastPatch`.
    """
    N = 4 * m * n
    edges: List[Tuple[int, int]] = []

    # intra-rhomb 4-cycle: 0-1-3-2-0
    for idx in range(m * n):
        b = idx * 4
        edges += [
            (b + 0, b + 1),
            (b + 1, b + 3),
            (b + 3, b + 2),
            (b + 2, b + 0),
        ]

    # inter-rhomb couplings (with optional PBC wrap)
    for r in range(m):
        row_base = r * n
        for c in range(n):
            b = (row_base + c) * 4
            if c < n - 1:  # right neighbour
                edges.append((b + 2, (row_base + c + 1) * 4 + 1))
            elif periodic and n > 1:
                edges.append((b + 2, row_base * 4 + 1))
            if r < m - 1:  # neighbour above
                edges.append((b + 3, ((r + 1) * n + c) * 4 + 0))
            elif periodic and m > 1:
                edges.append((b + 3, c * 4 + 0))

    pos = np.empty((N, 2), dtype=float)
    for r in range(m):
        for c in range(n):
            bx, by = c * _RHOMB_SPACING, r * _RHOMB_SPACING
            b = (r * n + c) * 4
            pos[b + 0] = (bx, by - _RHOMB_SIZE)
            pos[b + 1] = (bx - _RHOMB_SIZE, by)
            pos[b + 2] = (bx + _RHOMB_SIZE, by)
            pos[b + 3] = (bx, by + _RHOMB_SIZE)

    return N, edges, pos


# --- built-in hexagonal lattice (NetworkX-compatible) ------------------------
def _builtin_hex(
    m: int, n: int, periodic: bool = False
) -> Tuple[
    List[Coord], List[Tuple[Coord, Coord]], Dict[Coord, Tuple[float, float]]
]:
    """Reproduce ``networkx.hexagonal_lattice_graph(m, n, periodic)`` exactly.

    Returns coordinate-labelled ``(nodes, edges, pos)``.  Used as the scaffold
    for the kagome and tri-hexagonal lattices.  Periodic boundaries replicate
    NetworkX's ``contracted_nodes`` identifications via union-find.
    """
    if periodic and (n % 2 == 1 or m < 2 or n < 2):
        raise ValueError(
            "periodic hexagonal lattice needs m > 1, n > 1 and even n"
        )

    M = 2 * m
    nodes: set = set()
    edges: set = set()

    def add(a: Coord, b: Coord) -> None:
        edges.add((a, b))
        nodes.add(a)
        nodes.add(b)

    for i in range(n + 1):
        for j in range(M + 1):  # vertical (column) edges
            add((i, j), (i, j + 1))
    for i in range(n):
        for j in range(M + 2):  # horizontal (row) edges
            if i % 2 == j % 2:
                add((i, j), (i + 1, j))

    # drop the two single-edge corner nodes (NetworkX removes these)
    rm = {(0, M + 1), (n, (M + 1) * (n % 2))}
    nodes -= rm
    edges = {(a, b) for (a, b) in edges if a not in rm and b not in rm}

    if periodic:
        rep: Dict[Coord, Coord] = {nd: nd for nd in nodes}

        def find(x: Coord) -> Coord:
            while rep[x] != x:
                rep[x] = rep[rep[x]]
                x = rep[x]
            return x

        def merge(keep: Coord, drop: Coord) -> None:
            if keep in rep and drop in rep:
                rep[find(drop)] = find(keep)

        for i in range(n):  # (i, M)   -> (i, 0)
            merge((i, 0), (i, M))
        for i in range(1, n + 1):  # (i, M+1) -> (i, 1)
            merge((i, 1), (i, M + 1))
        for j in range(1, M):  # (n, j)   -> (0, j)
            merge((0, j), (n, j))

        merged_edges: set = set()
        for a, b in edges:
            ra, rb = find(a), find(b)
            if ra == rb:
                continue
            merged_edges.add((ra, rb) if ra <= rb else (rb, ra))
        nodes = {find(x) for x in nodes}
        edges = merged_edges

        last = (
            find((n, M)) if (n, M) in rep else None
        )  # NetworkX removes (n, M)
        if last is not None:
            nodes.discard(last)
            edges = {(a, b) for (a, b) in edges if a != last and b != last}

    def xpos(i: int, j: int) -> float:
        return 0.5 + i + i // 2 + (j % 2) * ((i % 2) - 0.5)

    pos = {(i, j): (xpos(i, j), _HEX_H * j) for (i, j) in nodes}
    return sorted(nodes), sorted(edges), pos


# --- kagome -------------------------------------------------------------------
def kagome(m: int, n: int, periodic: bool = False) -> GenResult:
    """Kagome lattice: line graph of the (built-in) hexagonal lattice.

    Bulk coordination ``z = 4``.  Each vertex is a hex *edge*; two kagome
    vertices are adjacent iff the underlying hex edges share an endpoint.
    """
    hex_nodes, hex_edges, hex_pos = _builtin_hex(m, n, periodic)

    # kagome vertex i  <->  hex edge hex_edges[i]
    edge_index = {e: i for i, e in enumerate(hex_edges)}
    N = len(hex_edges)

    # group hex edges by shared hex vertex, then clique each group
    incident: Dict[Coord, List[int]] = {nd: [] for nd in hex_nodes}
    for e, i in edge_index.items():
        incident[e[0]].append(i)
        incident[e[1]].append(i)

    edges: List[Tuple[int, int]] = []
    seen: set = set()
    for group in incident.values():
        for a in range(len(group)):
            for b in range(a + 1, len(group)):
                u, v = group[a], group[b]
                key = (u, v) if u < v else (v, u)
                if key not in seen:
                    seen.add(key)
                    edges.append(key)

    pos = np.empty((N, 2), dtype=float)
    for e, i in edge_index.items():
        (x0, y0), (x1, y1) = hex_pos[e[0]], hex_pos[e[1]]
        pos[i] = (0.5 * (x0 + x1), 0.5 * (y0 + y1))

    return N, edges, pos


# --- tri-hexagonal (3,12^2) ---------------------------------------------------
def tri_hexagonal(m: int, n: int, periodic: bool = False) -> GenResult:
    """Tri-hexagonal (truncated honeycomb) lattice; bulk coordination ``z = 3``.

    Each hex vertex ``v`` is truncated into a clique over its incident
    half-edges ``(v, nb)``; each hex edge ``(u, v)`` links ``(u, v)`` to
    ``(v, u)``.  Periodic boundaries are not supported (mirrors the NX code).
    """
    if periodic:
        raise NotImplementedError(
            "PBC not implemented for tri-hexagonal lattice"
        )

    hex_nodes, hex_edges, hex_pos = _builtin_hex(m, n, periodic=False)

    adj: Dict[Coord, List[Coord]] = {nd: [] for nd in hex_nodes}
    for a, b in hex_edges:
        adj[a].append(b)
        adj[b].append(a)

    # tri-hex vertex = directed half-edge (v, nb)
    half_edges: List[Tuple[Coord, Coord]] = []
    for v in hex_nodes:
        for nb in sorted(adj[v]):
            half_edges.append((v, nb))
    he_index = {he: i for i, he in enumerate(half_edges)}
    N = len(half_edges)

    edges: List[Tuple[int, int]] = []
    # triangle (clique) over each truncated vertex
    for v in hex_nodes:
        grp = [he_index[(v, nb)] for nb in sorted(adj[v])]
        for a in range(len(grp)):
            for b in range(a + 1, len(grp)):
                edges.append((grp[a], grp[b]))
    # bond connecting the two half-edges of every hex edge
    for a, b in hex_edges:
        edges.append((he_index[(a, b)], he_index[(b, a)]))

    pos = np.empty((N, 2), dtype=float)
    for (v, nb), i in he_index.items():
        xv, yv = hex_pos[v]
        xn, yn = hex_pos[nb]
        pos[i] = (
            xv + _TRIHEX_SHRINK * (xn - xv),
            yv + _TRIHEX_SHRINK * (yn - yv),
        )

    return N, edges, pos
