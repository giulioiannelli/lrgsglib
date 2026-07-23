from itertools import product
from typing import Tuple

from networkx import Graph

__all__ = [
    'generate_bcc_lattice',
    'generate_fcc_lattice'
]

# Canonical nearest-neighbour displacement sets (unit-cell coordinates),
# matching the GT engine bond-for-bond:
# BCC — every cell-centre atom bonds to the 8 surrounding corner atoms (z=8);
# FCC — every atom bonds to the 12 atoms half a face-diagonal away (z=12).
BCC_NN_DELTAS = [
    (dx, dy, dz)
    for dx in (0.5, -0.5) for dy in (0.5, -0.5) for dz in (0.5, -0.5)
]
FCC_NN_DELTAS = [
    (0.5, 0.5, 0), (0.5, -0.5, 0), (-0.5, 0.5, 0), (-0.5, -0.5, 0),
    (0.5, 0, 0.5), (0.5, 0, -0.5), (-0.5, 0, 0.5), (-0.5, 0, -0.5),
    (0, 0.5, 0.5), (0, 0.5, -0.5), (0, -0.5, 0.5), (0, -0.5, -0.5),
]


def _sublattice(dim, off, periodic: bool) -> list:
    """Nodes of one sublattice: unit-cell corners displaced by ``off``.

    Under open boundaries the sublattice loses its last unit cell along every
    axis its offset displaces (no partner cell to bond to) — the same
    convention as the GT engine, so both engines agree node-for-node.
    """
    ranges = [range(d if periodic or o == 0 else d - 1)
              for d, o in zip(dim, off)]
    return [(x + off[0], y + off[1], z + off[2])
            for x, y, z in product(*ranges)]


def _wrap(coord, dim):
    """Wrap a coordinate into the periodic box (exact for half-integers)."""
    return (coord[0] % dim[0], coord[1] % dim[1], coord[2] % dim[2])


def generate_bcc_lattice(dim, periodic: bool) -> Graph:
    """
    Generate a Body-Centered Cubic (BCC) lattice graph.

    Nearest-neighbour bonds only: every cell-centre atom connects to its 8
    surrounding corner atoms (coordination z=8 under PBC), identical to the
    GT engine's edge set.

    Parameters
    ----------
    dim : tuple of int
        Dimensions of the lattice (x, y, z).
    periodic : bool
        Whether to use periodic boundary conditions.

    Returns
    -------
    Graph
        A NetworkX graph representing the BCC lattice.
    """
    G = Graph()
    corners = _sublattice(dim, (0, 0, 0), periodic)
    centers = _sublattice(dim, (0.5, 0.5, 0.5), periodic)
    G.add_nodes_from(corners)
    G.add_nodes_from(centers)
    corner_set = set(corners)

    # every NN bond couples one centre and one corner, so generating from the
    # centre sublattice yields each edge exactly once
    edges = []
    for c in centers:
        for dx, dy, dz in BCC_NN_DELTAS:
            nb = (c[0] + dx, c[1] + dy, c[2] + dz)
            if periodic:
                nb = _wrap(nb, dim)
            if nb in corner_set:
                edges.append((c, nb, {'type': 'link'}))
    G.add_edges_from(edges)
    return G
#
def generate_fcc_lattice(dim: Tuple[int, int, int], periodic: bool) -> Graph:
    """
    Generate a Face-Centered Cubic (FCC) lattice graph.

    Nearest-neighbour bonds only: every atom connects to the 12 atoms half a
    face-diagonal away (coordination z=12 under PBC), identical to the GT
    engine's edge set.

    Parameters
    ----------
    dim : Tuple[int, int, int]
        Dimensions of the lattice (x, y, z).
    periodic : bool
        Whether to use periodic boundary conditions.

    Returns
    -------
    Graph
        A NetworkX graph representing the FCC lattice.
    """
    G = Graph()
    offsets = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    nodes = [n for off in offsets for n in _sublattice(dim, off, periodic)]
    G.add_nodes_from(nodes)
    node_set = set(nodes)

    # each undirected NN pair is generated from both endpoints; nx.Graph
    # deduplicates, mirroring the GT engine's added_edges set
    edges = []
    for u in nodes:
        for dx, dy, dz in FCC_NN_DELTAS:
            nb = (u[0] + dx, u[1] + dy, u[2] + dz)
            if periodic:
                nb = _wrap(nb, dim)
            if nb in node_set:
                edges.append((u, nb))
    G.add_edges_from(edges)
    return G
