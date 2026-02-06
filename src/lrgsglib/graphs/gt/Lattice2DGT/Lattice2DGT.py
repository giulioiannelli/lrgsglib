"""
Lattice2DGT: graph-tool-based 2D lattice signed graph.

This module re-exports the Lattice2DGT class from gt_patches for the
unified graph module structure.

The actual implementation remains in gt_patches.Lattice2DGT for backward
compatibility.

Example
-------
>>> from lrgsglib.graphs.gt.lattice import Lattice2DGT
>>> lat = Lattice2DGT(side1=50, geo='sqr', pflip=0.3)
>>> lat.flip_random_fract_edges()
>>> print(f"Nodes: {lat.N}, Edges: {lat.num_edges}")
"""

from ....gt_patches.Lattice2DGT import Lattice2DGT

__all__ = ["Lattice2DGT"]
