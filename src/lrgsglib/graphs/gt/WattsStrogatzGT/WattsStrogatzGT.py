"""
WattsStrogatzGT: graph-tool-based Watts-Strogatz small-world signed graph.

This module re-exports the WattsStrogatzGT class from gt_patches for the
unified graph module structure.

The actual implementation remains in gt_patches.random.watts_strogatz_gt for backward
compatibility.

Example
-------
>>> from lrgsglib.graphs.gt.random import WattsStrogatzGT
>>> ws = WattsStrogatzGT(n=100, k=4, p=0.3, pflip=0.1, seed=42)
>>> ws.flip_random_fract_edges()
>>> print(f"Nodes: {ws.N}, Edges: {ws.num_edges}")
"""

from ....gt_patches.random.watts_strogatz_gt import WattsStrogatzGT

__all__ = ["WattsStrogatzGT"]
