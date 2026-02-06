"""
Lattice3DGT: graph-tool-based 3D lattice signed graph.

This module re-exports the Lattice3DGT class from gt_patches for the
unified graph module structure.

The actual implementation remains in gt_patches.Lattice3DGT for backward
compatibility.

Example
-------
>>> from lrgsglib.graphs.gt.lattice import Lattice3DGT
>>> lat = Lattice3DGT(dim=10, geo='sc', pflip=0.3)
>>> lat.flip_random_fract_edges()
>>> print(f"Nodes: {lat.N}, Edges: {lat.num_edges}")
"""

from ....gt_patches.Lattice3DGT import Lattice3DGT

__all__ = ["Lattice3DGT"]
