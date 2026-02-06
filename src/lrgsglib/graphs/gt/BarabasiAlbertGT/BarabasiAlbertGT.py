"""
BarabasiAlbertGT: graph-tool-based Barabasi-Albert scale-free signed graph.

This module re-exports the BarabasiAlbertGT class from gt_patches for the
unified graph module structure.

The actual implementation remains in gt_patches.random.barabasi_albert_gt for backward
compatibility.

Example
-------
>>> from lrgsglib.graphs.gt.random import BarabasiAlbertGT
>>> ba = BarabasiAlbertGT(n=500, m=2, pflip=0.15, seed=42)
>>> ba.flip_random_fract_edges()
>>> print(f"Nodes: {ba.N}, Edges: {ba.num_edges}")
"""

from ....gt_patches.random.barabasi_albert_gt import BarabasiAlbertGT

__all__ = ["BarabasiAlbertGT"]
