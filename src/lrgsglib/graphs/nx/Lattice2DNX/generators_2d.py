from numpy import sqrt
from numpy.random import choice, rand
import numpy as np
from networkx import Graph, empty_graph, NetworkXError, set_node_attributes
from networkx.utils import pairwise
from networkx.utils.backends import _dispatchable as _nx_dispatchable, _registered_algorithms
from typing import Any
#
__all__ = [
    # "triangular_lattice_graph_modified",
    "triangular_lattice_graph_FastPatch",
    "hexagonal_lattice_graph_FastPatch",
    "squared_lattice_graph_FastPatch",
    "squared_lattice_SW_graph_FastPatch",
    "rhomb_octagonal_graph_FastPatch",
    "kagome_lattice_graph",
    "tri_hexagonal_lattice_graph",
]
#
def _dispatchable_safe(**kwargs):
    def decorator(func):
        name = kwargs.get("name", func.__name__)
        if name in _registered_algorithms:
            return func
        return _nx_dispatchable(func, **kwargs)
    return decorator


@_dispatchable_safe(graphs=None, returns_graph=True)
def triangular_lattice_graph_modified(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
    create_using: Any = None
) -> Graph:
    """
    Generates a triangular lattice graph with optional periodic boundary
    conditions (PBC). This function creates a graph representing a
    triangular lattice with `m` rows and `n` columns. Nodes in the lattice
    are connected in a manner that forms a pattern of equilateral triangles.
    When periodic boundary conditions are enabled, the lattice simulates a
    toroidal surface where edges wrap around the opposite side, creating a
    continuous pattern.

    Parameters
    ----------
    m : int
        The number of rows in the lattice. Must be a positive integer.
    n : int
        The number of columns in the lattice. Must be a positive integer.
    periodic : bool, optional (default=False)
        If True, applies periodic boundary conditions, simulating a toroidal
        topology. Requires `m >= 3` and `n >= 5`.
    with_positions : bool, optional (default=True)
        If True, calculates and stores the positions of each node in the
        node attribute 'pos', arranging nodes in equilateral triangles.
    create_using : Any, optional
        A NetworkX graph constructor. If None, a default graph is created.

    Returns
    -------
    Graph
        A NetworkX graph object representing the m by n triangular lattice,
        optionally with PBC.

    Raises
    ------
    ValueError
        If `m` or `n` is not a positive integer.
    NetworkXError
        If periodic is True and `m < 3` or `n < 5`.
    """
    if not isinstance(m, int) or m <= 0:
        raise ValueError("`m` must be a positive integer.")
    if not isinstance(n, int) or n <= 0:
        raise ValueError("`n` must be a positive integer.")

    G = empty_graph(0, create_using)
    if n == 0 or m == 0:
        return G

    if periodic:
        if n < 5 or m < 3:
            msg = f"m > 2 and n > 4 required for periodic. m={m}, n={n}"
            raise NetworkXError(msg)

    rows = range(m)
    cols = range(n)

    # Identify boundary nodes if periodic
    if periodic:
        for i in range(n + 1):
            for j in range(m + 1):
                G.add_node((i % n, j % m))  # Add node with PBC
                # Add horizontal edges within the grid, with PBC for the last
                # column
                if i < n or j % 2 == 0:  # For even rows, wrap horizontally
                    G.add_edge((i % n, j % m), ((i + 1) % n, j % m))
                # Add vertical and diagonal edges, with PBC for the last row
                if j < m:
                    G.add_edge((i % n, j % m), (i % n, (j + 1) % m))
                    if j % 2:  # Diagonal for even rows
                        G.add_edge(
                            (i % n, j % m), ((i + 1) % n, (j + 1) % m)
                        )
                    else:  # Diagonal for odd rows, wrapping if at the edge
                        G.add_edge(
                            (i % n, j % m), ((i - 1) % n, (j + 1) % m)
                        )
    else:
        # Make grid
        G.add_edges_from(
            ((i, j), (i + 1, j)) for j in rows for i in cols[: n - 1]
        )
        G.add_edges_from(
            ((i, j), (i, j + 1)) for j in rows[: m - 1] for i in cols
        )
        # Add diagonals
        G.add_edges_from(
            ((i, j), (i + 1, j + 1))
            for j in rows[1 : m - 1 : 2]
            for i in cols[: n - 1]
        )
        G.add_edges_from(
            ((i + 1, j), (i, j + 1))
            for j in rows[: m - 1 : 2]
            for i in cols[: n - 1]
        )

    # Add position node attributes
    if with_positions:
        ii = (i for i in cols for j in rows)
        jj = (j for i in cols for j in rows)
        xx = (0.5 * (j % 2) + i for i in cols for j in rows)
        h = sqrt(3) / 2
        if periodic:
            yy = (h * j + 0.01 * i * i for i in cols for j in rows)
        else:
            yy = (h * j for i in cols for j in rows)
        pos = {
            (i, j): (x, y)
            for i, j, x, y in zip(ii, jj, xx, yy)
            if (i, j) in G
        }
        set_node_attributes(G, pos, "pos")

    return G
#
@_dispatchable_safe(graphs=None, returns_graph=True)
def triangular_lattice_graph_FastPatch(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
    create_using: Any = None,
    bend_positions: bool = False
) -> Graph:
    """
    Generates a triangular lattice graph with optional periodic boundary
    conditions (PBC). This function creates a graph representing a
    triangular lattice with `m` rows and `n` columns. Nodes in the lattice
    are connected in a manner that forms a pattern of equilateral triangles.
    When periodic boundary conditions are enabled, the lattice simulates a
    toroidal surface where edges wrap around the opposite side, creating a
    continuous pattern.

    Parameters
    ----------
    m : int
        The number of rows in the lattice. Must be a positive integer.
    n : int
        The number of columns in the lattice. Must be a positive integer.
    periodic : bool, optional (default=False)
        If True, applies periodic boundary conditions, simulating a toroidal
        topology. Requires `m >= 3` and `n >= 5`.
    with_positions : bool, optional (default=True)
        If True, calculates and stores the positions of each node in the
        node attribute 'pos', arranging nodes in equilateral triangles.
    create_using : Any, optional
        A NetworkX graph constructor. If None, a default graph is created.
    bend_positions : bool, optional (default=False)
        If True, applies a slight bend to the positions for visualization.

    Returns
    -------
    Graph
        A NetworkX graph object representing the m by n triangular lattice,
        optionally with PBC.

    Raises
    ------
    ValueError
        If `m` or `n` is not a positive integer.
    NetworkXError
        If periodic is True and `m < 3` or `n < 5`.
    """

    if not isinstance(m, int) or m <= 0:
        raise ValueError("`m` must be a positive integer.")
    if not isinstance(n, int) or n <= 0:
        raise ValueError("`n` must be a positive integer.")

    G = empty_graph(0, create_using)

    if periodic:
        if n < 5 or m < 3:
            msg = f"m > 2 and n > 4 required for periodic. m={m}, n={n}"
            raise NetworkXError(msg)

    rows = range(m)
    cols = range(n)

    # Identify boundary nodes if periodic
    if periodic:
        for i in range(n + 1):
            for j in range(m + 1):
                G.add_node((i % n, j % m))  # Add node with PBC
                # Add horizontal edges within the grid, with PBC for the last
                # column
                if i < n or j % 2 == 0:  # For even rows, wrap horizontally
                    G.add_edge((i % n, j % m), ((i + 1) % n, j % m))
                # Add vertical and diagonal edges, with PBC for the last row
                if j < m:
                    G.add_edge((i % n, j % m), (i % n, (j + 1) % m))
                    if j % 2:  # Diagonal for even rows
                        G.add_edge(
                            (i % n, j % m), ((i + 1) % n, (j + 1) % m)
                        )
                    else:  # Diagonal for odd rows, wrapping if at the edge
                        G.add_edge(
                            (i % n, j % m), ((i - 1) % n, (j + 1) % m)
                        )
    else:
        # Make grid
        G.add_edges_from(
            ((i, j), (i + 1, j)) for j in rows for i in cols[: n - 1]
        )
        G.add_edges_from(
            ((i, j), (i, j + 1)) for j in rows[: m - 1] for i in cols
        )
        # Add diagonals
        G.add_edges_from(
            ((i, j), (i + 1, j + 1))
            for j in rows[1 : m - 1 : 2]
            for i in cols[: n - 1]
        )
        G.add_edges_from(
            ((i + 1, j), (i, j + 1))
            for j in rows[: m - 1 : 2]
            for i in cols[: n - 1]
        )

    # Add position node attributes
    if with_positions:
        ii = (i for i in rows for j in cols)
        jj = (j for i in rows for j in cols)
        xx = (0.5 * (i % 2) + j for i in rows for j in cols)
        h = sqrt(3) / 2
        if periodic and bend_positions:
            yy = (h * i + 0.01 * j * j for i in rows for j in cols)
        else:
            yy = (h * i for i in rows for j in cols)
        pos = {
            (j, i): (x, y)
            for i, j, x, y in zip(ii, jj, xx, yy)
            if (j, i) in G
        }
        set_node_attributes(G, pos, "pos")

    return G
#
@_dispatchable_safe(graphs=None, returns_graph=True)
def hexagonal_lattice_graph_FastPatch(
    n: int,
    m: int,
    periodic: bool = False,
    with_positions: bool = True,
    create_using: Any = None,
    bend_positions: bool = False
) -> Graph:
    """
    Generate a hexagonal lattice graph with optional periodic boundary
    conditions.

    Parameters
    ----------
    n : int
        Number of columns of hexagons. Must be a positive integer.
    m : int
        Number of rows of hexagons. Must be a positive integer.
    periodic : bool, optional (default=False)
        If True, applies periodic boundary conditions to create a toroidal
        topology.
    with_positions : bool, optional (default=True)
        If True, calculates and assigns positions to nodes to represent a
        hexagonal lattice.
    create_using : Any, optional
        A NetworkX graph constructor. If None, a default graph is created.
    bend_positions : bool, optional (default=False)
        If True, applies a slight bend to the positions for visualization.

    Returns
    -------
    Graph
        A NetworkX graph object representing the hexagonal lattice.

    Raises
    ------
    ValueError
        If `n` or `m` is not a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("`n` must be a positive integer.")
    if not isinstance(m, int) or m <= 0:
        raise ValueError("`m` must be a positive integer.")

    G = empty_graph(0, create_using)

    rows = range(m)
    cols = range(n)

    # Add edges for the hexagonal lattice
    col_edges = (((i, j), (i, j + 1)) for i in cols for j in rows[: m - 1])
    row_edges = (
        ((i, j), (i + 1, j))
        for i in cols[: n - 1]
        for j in rows
        if i % 2 == j % 2
    )
    G.add_edges_from(col_edges)
    G.add_edges_from(row_edges)

    if periodic:
        more_row_edges = (
            ((cols[0], j), (cols[n - 1], j)) for j in rows if j % 2
        )
        more_col_edges = (
            ((i, rows[0]), (i, rows[m - 1])) for i in cols if (i + 1) % 2
        )
        G.add_edges_from(more_row_edges)
        G.add_edges_from(more_col_edges)

    # Optionally assign positions for visualization
    if with_positions:
        ii = (i for i in cols for j in rows)
        jj = (j for i in cols for j in rows)
        xx = (
            0.5 + i + i // 2 + (j % 2) * ((i % 2) - 0.5)
            for i in cols
            for j in rows
        )
        h = sqrt(3) / 2
        if periodic and bend_positions:
            yy = (h * j + 0.01 * i * i for i in cols for j in rows)
        else:
            yy = (h * j for i in cols for j in rows)
        pos = {
            (i, j): (x, y)
            for i, j, x, y in zip(ii, jj, xx, yy)
            if (i, j) in G
        }
        set_node_attributes(G, pos, "pos")

    return G
#
@_dispatchable_safe(graphs=None, returns_graph=True)
def squared_lattice_graph_FastPatch(
    m: int,
    n: int,
    periodic: bool = False,
    create_using: Any = None,
    with_positions: bool = True,
    bend_positions: bool = False
) -> Graph:
    """
    Returns the two-dimensional grid graph.

    The grid graph has each node connected to its four nearest neighbors.

    Parameters
    ----------
    m : int
        Number of rows in the grid. Must be a positive integer.
    n : int
        Number of columns in the grid. Must be a positive integer.
    periodic : bool or iterable, optional (default=False)
        If `periodic` is True, both dimensions are periodic. If False, none
        are periodic. If `periodic` is iterable, it should yield 2 bool
        values indicating whether the 1st and 2nd axes, respectively, are
        periodic.
    create_using : Any, optional
        A NetworkX graph constructor. If None, a default graph is created.
    with_positions : bool, optional (default=True)
        If True, calculates and assigns positions to nodes for visualization.
    bend_positions : bool, optional (default=False)
        If True, applies a slight bend to the positions for visualization.

    Returns
    -------
    Graph
        The (possibly periodic) grid graph of the specified dimensions.

    Raises
    ------
    ValueError
        If `m` or `n` is not a positive integer.
    """
    if not isinstance(m, int) or m <= 0:
        raise ValueError("`m` must be a positive integer.")
    if not isinstance(n, int) or n <= 0:
        raise ValueError("`n` must be a positive integer.")

    G = empty_graph(0, create_using)
    rows = range(m)
    cols = range(n)
    G.add_nodes_from((i, j) for i in rows for j in cols)
    G.add_edges_from(((i, j), (pi, j)) for pi, i in pairwise(rows) 
                     for j in cols)
    G.add_edges_from(((i, j), (i, pj)) for i in rows 
                     for pj, j in pairwise(cols))

    try:
        periodic_r, periodic_c = periodic
    except TypeError:
        periodic_r = periodic_c = periodic

    if periodic_r and len(rows) > 2:
        first = rows[0]
        last = rows[-1]
        G.add_edges_from(((first, j), (last, j)) for j in cols)
    if periodic_c and len(cols) > 2:
        first = cols[0]
        last = cols[-1]
        G.add_edges_from(((i, first), (i, last)) for i in rows)

    if G.is_directed():
        G.add_edges_from((v, u) for u, v in G.edges())

    if with_positions:
        ii = (i for i in rows for j in cols)
        jj = (j for i in rows for j in cols)
        xx = (i for i in rows for j in cols)
        yy = (j for i in rows for j in cols)
        if periodic and bend_positions:
            xx = (i + 0.02 * j * j for i in rows for j in cols)
            yy = (j + 0.02 * i * i for i in rows for j in cols)
        pos = {
            (i, j): (x, y)
            for i, j, x, y in zip(ii, jj, xx, yy)
            if (i, j) in G
        }
        set_node_attributes(G, pos, "pos")

    return G
#
@_dispatchable_safe(graphs=None, returns_graph=True)
def squared_lattice_SW_graph_FastPatch(
    m: int,
    n: int,
    prew: float = 0,
    periodic: bool = False,
    create_using: Any = None,
    with_positions: bool = True
) -> Graph:
    """
    Generates a small-world squared lattice graph with optional periodic
    boundary conditions and rewiring of edges.

    Parameters
    ----------
    m : int
        Number of rows in the grid. Must be a positive integer.
    n : int
        Number of columns in the grid. Must be a positive integer.
    prew : float, optional (default=0)
        Probability of rewiring each edge. Must be in the range [0, 1].
    periodic : bool, optional (default=False)
        If True, applies periodic boundary conditions to the grid.
    create_using : Any, optional
        A NetworkX graph constructor. If None, a default graph is created.
    with_positions : bool, optional (default=True)
        If True, calculates and assigns positions to nodes for visualization.

    Returns
    -------
    Graph
        A NetworkX graph object representing the small-world squared lattice.

    Raises
    ------
    ValueError
        If `m` or `n` is not a positive integer, or if `prew` is not in
        the range [0, 1].
    """
    if not isinstance(m, int) or m <= 0:
        raise ValueError("`m` must be a positive integer.")
    if not isinstance(n, int) or n <= 0:
        raise ValueError("`n` must be a positive integer.")
    if not (0 <= prew <= 1):
        raise ValueError("`prew` must be in the range [0, 1].")

    G = squared_lattice_graph_FastPatch(
        m, n, periodic, create_using, with_positions
    )

    nodes = list(G.nodes())
    edges = list(G.edges())
    num_edges = len(edges)

    rewire_flags = rand(num_edges) < prew
    new_neighbors = choice(len(nodes), size=num_edges, replace=True)

    for i, edge in enumerate(edges):
        if rewire_flags[i]:
            u, v = edge
            fixed_node = u
            G.remove_edge(u, v)

            new_neigh = nodes[new_neighbors[i]]
            while new_neigh == fixed_node or G.has_edge(fixed_node, new_neigh):
                new_neigh = nodes[choice(len(nodes))]

            G.add_edge(fixed_node, new_neigh)

    return G

@_dispatchable_safe(graphs=None, returns_graph=True)
def rhomb_octagonal_graph_FastPatch(
    m: int,
    n: int,
    periodic: bool = False,
    create_using: Any = None,
    with_positions: bool = True,
    bend_positions: bool = False
) -> Graph:
    """
    Generates a rhomb-octagonal lattice graph where each node of a square 
    lattice is replaced by four nodes forming a rhomb. This creates a 
    dual-scale lattice with rhombs at one scale and octagonal cells at 
    another scale.

    Parameters
    ----------
    m : int
        Number of rows of rhombs in the lattice. Must be a positive integer.
    n : int
        Number of columns of rhombs in the lattice. Must be a positive integer.
        The total number of nodes will be 4*m*n.
    periodic : bool, optional (default=False)
        If True, applies periodic boundary conditions to the lattice.
    create_using : Any, optional
        A NetworkX graph constructor. If None, a default graph is created.
    with_positions : bool, optional (default=True)
        If True, calculates and assigns positions to nodes for visualization.
    bend_positions : bool, optional (default=False)
        If True, applies a slight bend to the positions for visualization.

    Returns
    -------
    Graph
        A NetworkX graph object representing the rhomb-octagonal lattice.

    Raises
    ------
    ValueError
        If `m` or `n` is not a positive integer.
    """
    if not isinstance(m, int) or m <= 0:
        raise ValueError("`m` must be a positive integer.")
    if not isinstance(n, int) or n <= 0:
        raise ValueError("`n` must be a positive integer.")

    N = 4 * m * n
    G = empty_graph(0, create_using)

    # Add all nodes
    G.add_nodes_from(range(N))

    # Add edges within each rhomb
    # Each rhomb has nodes arranged as: 3(top), 1(right), 0(bottom), 2(left)
    # Connect each node to its two neighbors in the rhomb cycle: 0-1-3-2-0
    intra_edges = []
    for rhomb_idx in range(m * n):
        base = rhomb_idx * 4
        # Intra-rhomb edges forming a 4-cycle: 0-1-3-2-0
        intra_edges.extend(
            [
                (base + 0, base + 1),  # bottom to right
                (base + 1, base + 3),  # right to top
                (base + 3, base + 2),  # top to left
                (base + 2, base + 0),  # left to bottom
            ]
        )
    G.add_edges_from(intra_edges)

    # Add inter-rhomb edges
    inter_edges = []
    for rhomb_row in range(m):
        row_base = rhomb_row * n
        for rhomb_col in range(n):
            base = (row_base + rhomb_col) * 4

            # Connect to rhomb to the right (horizontally)
            if rhomb_col < n - 1:
                right_rhomb_base = (row_base + rhomb_col + 1) * 4
                inter_edges.append(
                    (base + 2, right_rhomb_base + 1)
                )  # right to left
            elif periodic and n > 1:
                leftmost_rhomb_base = row_base * 4
                inter_edges.append((base + 2, leftmost_rhomb_base + 1))

            # Connect to rhomb above (vertically)
            if rhomb_row < m - 1:
                above_rhomb_base = ((rhomb_row + 1) * n + rhomb_col) * 4
                inter_edges.append((base + 3, above_rhomb_base + 0))
            elif periodic and m > 1:
                bottom_rhomb_base = rhomb_col * 4
                inter_edges.append((base + 3, bottom_rhomb_base + 0))
    G.add_edges_from(inter_edges)

    # Add position attributes if requested
    if with_positions:
        pos = {}

        # Scale factor for spacing between rhombs
        rhomb_size = 0.4  # Half-distance from center to edge
        rhomb_spacing = 4.0 * rhomb_size  # Rhombs touch edge-to-edge

        for rhomb_row in range(m):
            for rhomb_col in range(n):
                base_x = rhomb_col * rhomb_spacing
                base_y = rhomb_row * rhomb_spacing

                if periodic and bend_positions:
                    base_x += 0.02 * rhomb_row * rhomb_row
                    base_y += 0.02 * rhomb_col * rhomb_col

                base = (rhomb_row * n + rhomb_col) * 4
                pos[base + 0] = (base_x, base_y - rhomb_size)
                pos[base + 1] = (base_x - rhomb_size, base_y)
                pos[base + 2] = (base_x + rhomb_size, base_y)
                pos[base + 3] = (base_x, base_y + rhomb_size)

        set_node_attributes(G, pos, "pos")

    return G


# ---- Kagome & tri-hexagonal generators ----

import networkx as _nx


def kagome_lattice_graph(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
    bend_positions: bool = False,
) -> Graph:
    """Generate a kagome lattice (line graph of the hexagonal lattice).

    The kagome lattice has coordination z=4 in the bulk.
    It appears in Fujii & Tokunaga (2012) with BPT(primal)=0.524.

    Parameters
    ----------
    m, n : int
        Side parameters for the underlying hexagonal lattice.
    periodic : bool
        Periodic boundary conditions.
    with_positions : bool
        Store node positions.
    bend_positions : bool
        Ignored (kept for API compatibility).

    Returns
    -------
    Graph
        Kagome lattice with coordinate-tuple node labels.
    """
    hex_graph = _nx.hexagonal_lattice_graph(m, n, periodic=periodic)
    G = _nx.line_graph(hex_graph)

    if with_positions:
        hex_pos = _nx.get_node_attributes(hex_graph, "pos")
        pos = {}
        for node in G.nodes():
            u, v = node
            if u in hex_pos and v in hex_pos:
                xu, yu = hex_pos[u]
                xv, yv = hex_pos[v]
                pos[node] = ((xu + xv) / 2.0, (yu + yv) / 2.0)
        set_node_attributes(G, pos, "pos")

    return G


def tri_hexagonal_lattice_graph(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
    bend_positions: bool = False,
) -> Graph:
    r"""Generate a (3,12^2) tri-hexagonal lattice.

    Constructed by truncating each vertex of the hexagonal lattice:
    every vertex becomes a triangle, every hexagonal face becomes a
    dodecagon.  Bulk coordination is z=3.  Appears in Fujii &
    Tokunaga (2012) with BPT(primal)=0.740.

    Parameters
    ----------
    m, n : int
        Side parameters for the underlying hexagonal lattice.
    periodic : bool
        Periodic boundary conditions (not yet supported).
    with_positions : bool
        Store node positions.
    bend_positions : bool
        Ignored.

    Returns
    -------
    Graph
        Tri-hexagonal lattice with coordinate-tuple node labels.
    """
    if periodic:
        raise NotImplementedError(
            "PBC not yet implemented for tri-hexagonal lattice"
        )

    hex_graph = _nx.hexagonal_lattice_graph(m, n, periodic=False)
    G = Graph()

    for v in hex_graph.nodes():
        nbrs = sorted(hex_graph.neighbors(v))
        tri_nodes = [(v, nb) for nb in nbrs]
        for tn in tri_nodes:
            G.add_node(tn)
        for i in range(len(tri_nodes)):
            for j in range(i + 1, len(tri_nodes)):
                G.add_edge(tri_nodes[i], tri_nodes[j])

    for u, v in hex_graph.edges():
        G.add_edge((u, v), (v, u))

    if with_positions:
        hex_pos = _nx.get_node_attributes(hex_graph, "pos")
        pos = {}
        shrink = 0.25
        for v in hex_graph.nodes():
            nbrs = sorted(hex_graph.neighbors(v))
            if v not in hex_pos:
                continue
            xv, yv = hex_pos[v]
            for nb in nbrs:
                if nb not in hex_pos:
                    continue
                xn, yn = hex_pos[nb]
                pos[(v, nb)] = (
                    xv + shrink * (xn - xv),
                    yv + shrink * (yn - yv),
                )
        set_node_attributes(G, pos, "pos")

    return G
