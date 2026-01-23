"""
Consolidated lattice graph generators.

This module provides a unified interface to all lattice generation functions
used across the nx_patches module, supporting 2D, 3D, and ND lattices.

2D Generators
-------------
- squared_lattice_graph: Square lattice with optional PBC
- triangular_lattice_graph: Triangular lattice with optional PBC
- hexagonal_lattice_graph: Hexagonal (honeycomb) lattice with optional PBC
- squared_lattice_SW_graph: Small-world square lattice
- rhomb_octagonal_graph: Rhomb-octagonal lattice

3D Generators
-------------
- simple_cubic_lattice: Simple cubic lattice
- bcc_lattice: Body-centered cubic lattice
- fcc_lattice: Face-centered cubic lattice

Example
-------
>>> from lrgsglib.nx_patches._generators import lattice
>>> G = lattice.squared_lattice_graph(10, 10, periodic=True)
>>> G.number_of_nodes()
100
"""

from networkx import Graph, grid_2d_graph

# Import 2D generators from existing module
from ..Lattice2D.generators_2d import (
    triangular_lattice_graph_FastPatch,
    hexagonal_lattice_graph_FastPatch,
    squared_lattice_graph_FastPatch,
    squared_lattice_SW_graph_FastPatch,
    rhomb_octagonal_graph_FastPatch,
)

# Import 3D generators from existing module
from ..Lattice3D.generators_3d import (
    generate_bcc_lattice,
    generate_fcc_lattice,
)

__all__ = [
    # 2D lattice generators
    "squared_lattice_graph",
    "triangular_lattice_graph",
    "hexagonal_lattice_graph",
    "squared_lattice_SW_graph",
    "rhomb_octagonal_graph",
    # 3D lattice generators
    "simple_cubic_lattice",
    "bcc_lattice",
    "fcc_lattice",
    # Legacy names (for backward compatibility)
    "squared_lattice_graph_FastPatch",
    "triangular_lattice_graph_FastPatch",
    "hexagonal_lattice_graph_FastPatch",
    "squared_lattice_SW_graph_FastPatch",
    "rhomb_octagonal_graph_FastPatch",
    "generate_bcc_lattice",
    "generate_fcc_lattice",
]


# =============================================================================
# 2D Lattice Generators (clean interface)
# =============================================================================

def squared_lattice_graph(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
) -> Graph:
    """
    Generate a square lattice graph.

    Parameters
    ----------
    m : int
        Number of rows.
    n : int
        Number of columns.
    periodic : bool, default False
        Whether to use periodic boundary conditions.
    with_positions : bool, default True
        Whether to compute node positions.

    Returns
    -------
    Graph
        A NetworkX graph representing the square lattice.

    Example
    -------
    >>> G = squared_lattice_graph(5, 5, periodic=True)
    >>> G.number_of_nodes()
    25
    """
    return squared_lattice_graph_FastPatch(m, n, periodic=periodic, with_positions=with_positions)


def triangular_lattice_graph(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
) -> Graph:
    """
    Generate a triangular lattice graph.

    Parameters
    ----------
    m : int
        Number of rows.
    n : int
        Number of columns.
    periodic : bool, default False
        Whether to use periodic boundary conditions.
    with_positions : bool, default True
        Whether to compute node positions.

    Returns
    -------
    Graph
        A NetworkX graph representing the triangular lattice.
    """
    return triangular_lattice_graph_FastPatch(m, n, periodic=periodic, with_positions=with_positions)


def hexagonal_lattice_graph(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
) -> Graph:
    """
    Generate a hexagonal (honeycomb) lattice graph.

    Parameters
    ----------
    m : int
        Number of rows.
    n : int
        Number of columns.
    periodic : bool, default False
        Whether to use periodic boundary conditions.
    with_positions : bool, default True
        Whether to compute node positions.

    Returns
    -------
    Graph
        A NetworkX graph representing the hexagonal lattice.
    """
    return hexagonal_lattice_graph_FastPatch(m, n, periodic=periodic, with_positions=with_positions)


def squared_lattice_SW_graph(
    m: int,
    n: int,
    p: float = 0.0,
    periodic: bool = False,
    with_positions: bool = True,
) -> Graph:
    """
    Generate a small-world square lattice graph.

    Parameters
    ----------
    m : int
        Number of rows.
    n : int
        Number of columns.
    p : float, default 0.0
        Probability of rewiring each edge.
    periodic : bool, default False
        Whether to use periodic boundary conditions.
    with_positions : bool, default True
        Whether to compute node positions.

    Returns
    -------
    Graph
        A NetworkX graph representing the small-world square lattice.
    """
    return squared_lattice_SW_graph_FastPatch(m, n, p=p, periodic=periodic, with_positions=with_positions)


def rhomb_octagonal_graph(
    m: int,
    n: int,
    periodic: bool = False,
    with_positions: bool = True,
) -> Graph:
    """
    Generate a rhomb-octagonal lattice graph.

    Parameters
    ----------
    m : int
        Number of rows.
    n : int
        Number of columns.
    periodic : bool, default False
        Whether to use periodic boundary conditions.
    with_positions : bool, default True
        Whether to compute node positions.

    Returns
    -------
    Graph
        A NetworkX graph representing the rhomb-octagonal lattice.
    """
    return rhomb_octagonal_graph_FastPatch(m, n, periodic=periodic, with_positions=with_positions)


# =============================================================================
# 3D Lattice Generators (clean interface)
# =============================================================================

def simple_cubic_lattice(
    dim: tuple[int, int, int],
    periodic: bool = False,
) -> Graph:
    """
    Generate a simple cubic lattice graph.

    Parameters
    ----------
    dim : tuple[int, int, int]
        Dimensions (x, y, z) of the lattice.
    periodic : bool, default False
        Whether to use periodic boundary conditions.

    Returns
    -------
    Graph
        A NetworkX graph representing the simple cubic lattice.

    Example
    -------
    >>> G = simple_cubic_lattice((3, 3, 3), periodic=True)
    >>> G.number_of_nodes()
    27
    """
    return grid_2d_graph(*dim[:2]) if len(dim) == 2 else _simple_cubic_impl(dim, periodic)


def _simple_cubic_impl(dim: tuple[int, int, int], periodic: bool) -> Graph:
    """Internal implementation for simple cubic lattice."""
    from networkx import grid_graph
    if periodic:
        return grid_graph(dim=dim, periodic=True)
    return grid_graph(dim=dim, periodic=False)


def bcc_lattice(
    dim: tuple[int, int, int],
    periodic: bool = False,
) -> Graph:
    """
    Generate a body-centered cubic (BCC) lattice graph.

    Parameters
    ----------
    dim : tuple[int, int, int]
        Dimensions (x, y, z) of the lattice.
    periodic : bool, default False
        Whether to use periodic boundary conditions.

    Returns
    -------
    Graph
        A NetworkX graph representing the BCC lattice.
    """
    return generate_bcc_lattice(dim, periodic)


def fcc_lattice(
    dim: tuple[int, int, int],
    periodic: bool = False,
) -> Graph:
    """
    Generate a face-centered cubic (FCC) lattice graph.

    Parameters
    ----------
    dim : tuple[int, int, int]
        Dimensions (x, y, z) of the lattice.
    periodic : bool, default False
        Whether to use periodic boundary conditions.

    Returns
    -------
    Graph
        A NetworkX graph representing the FCC lattice.
    """
    return generate_fcc_lattice(dim, periodic)
