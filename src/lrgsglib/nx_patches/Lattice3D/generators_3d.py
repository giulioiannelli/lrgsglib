from networkx import Graph
from typing import Tuple
#
from ...utils.basic.iterables import cProd_Iter, cProd_Iter_adj, cProdSel_Iter

__all__ = [
    'generate_bcc_lattice',
    'generate_fcc_lattice'
]

def generate_bcc_lattice(dim, periodic: bool) -> Graph:
    """
    Generate a Body-Centered Cubic (BCC) lattice graph.

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
    offsets = [(0, 0, 0), (0.5, 0.5, 0.5)]
    offsets2 = [(0.5, 0.5, 0.5), (0.5, 0.5, -0.5),
                (0.5, -0.5, 0.5), (-0.5, 0.5, 0.5),
                (0.5, -0.5, -0.5), (-0.5, -0.5, -0.5),
                (-0.5, 0.5, -0.5), (-0.5, -0.5, 0.5)]
    range_adjust = 0 if periodic else -1

    nodes = [(x + ox, y + oy, z + oz) 
             for ox, oy, oz in offsets
             for x, y, z in cProd_Iter_adj(dim, range_adjust)]
    G.add_nodes_from(nodes)

    edges = [((x + dx, y + dy, z + dz), 
              (x + ddx + dx, y + ddy + dy, z + ddz + dz), 
              {'type': 'link'})
             for x, y, z in cProd_Iter(dim) 
             for dx, dy, dz in offsets
             for ddx, ddy, ddz in offsets2
             if (x + ddx + dx, y + ddy + dy, z + ddz + dz) in G.nodes()]

    G.add_edges_from(edges)
    edges = [((x + dx, y + dy, z + dz), 
              (x + ddx + dx, y + ddy + dy, z + ddz + dz), 
              {'type': 'box'})
             for x, y, z in cProd_Iter(dim) 
             for dx, dy, dz in offsets
             for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
             if (x + ddx + dx, y + ddy + dy, z + ddz + dz) in G.nodes()]

    G.add_edges_from(edges)

    if periodic:
        # For BCC lattice PBC
        G.add_edges_from([((x, y, 0), 
                          (x+0.5, y+0.5, dim[2] - 0.5), 
                          {'type': 'link'}) 
                          for x, y in cProdSel_Iter(dim, (0, 1))])
        G.add_edges_from([((x, 0, 0), 
                          (x+0.5, dim[1]-0.5, dim[2] - 0.5), 
                          {'type': 'link'}) 
                          for x in range(dim[0])])
        G.add_edges_from([((0, y, 0), 
                          (dim[0]-0.5, y+0.5, dim[2] - 0.5), 
                          {'type': 'link'}) 
                          for y in range(dim[1])])
        G.add_edges_from([((0, 0, z), 
                          (dim[0]-0.5, dim[1]-0.5, z+0.5), 
                          {'type': 'link'}) 
                          for z in range(dim[2])])
        G.add_edges_from([((x, 0, z), 
                          (x+0.5, dim[1]-0.5, z+0.5), 
                          {'type': 'link'}) 
                          for x, z in cProdSel_Iter(dim, (0, 2))])
        G.add_edges_from([((0, y, z), 
                          (dim[0]-0.5, y+0.5, z+0.5), 
                          {'type': 'link'}) 
                          for y, z in cProdSel_Iter(dim, (1, 2))])

    return G
#
def generate_fcc_lattice(dim: Tuple[int, int, int], periodic: bool) -> Graph:
    """
    Generate a Face-Centered Cubic (FCC) lattice graph.

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
    range_adjust = 0 if periodic else -1

    nodes = [(x + ox, y + oy, z + oz) 
             for ox, oy, oz in offsets
             for x, y, z in cProd_Iter_adj(dim, range_adjust)]
    G.add_nodes_from(nodes)

    edges = [((x + dx, y + dy, z + dz), 
              (x + ddx + dx, y + ddy + dy, z + ddz + dz))
             for x, y, z in cProd_Iter(dim) 
             for dx, dy, dz in offsets
             for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
             if (x + ddx + dx, y + ddy + dy, z + ddz + dz) in G.nodes()]

    G.add_edges_from(edges)

    if periodic:
        # For FCC lattice periodic boundary conditions
        G.add_edges_from([((x + ox, y + oy, oz), 
                          (x + ox, y + oy, dim[2] - 1 + oz)) 
                          for ox, oy, oz in offsets
                          for x, y in cProdSel_Iter(dim, (0, 1))])
        G.add_edges_from([((x + ox, oy, z + oz), 
                          (x + ox, dim[1] - 1 + oy, z + oz))
                          for ox, oy, oz in offsets
                          for x, z in cProdSel_Iter(dim, (0, 2))])
        G.add_edges_from([((ox, y + oy, z + oz), 
                          (dim[0] - 1 + ox, y + oy, z + oz))
                          for ox, oy, oz in offsets
                          for y, z in cProdSel_Iter(dim, (1, 2))])

    return G