"""
SignedGraphGT: graph-tool-based signed graph implementation.

This module re-exports the SignedGraphGT class from gt_patches for the
unified graph module structure.

The actual implementation remains in gt_patches.signed_graph_gt for backward
compatibility.

Example
-------
>>> from lrgsglib.graphs.gt.base import SignedGraphGT
>>> sg = SignedGraphGT.from_complete(n=50, pflip=0.2, seed=42)
>>> sg.flip_random_fract_edges()
>>> print(f"Nodes: {sg.N}, Edges: {sg.num_edges}")
"""

from ....gt_patches.signed_graph_gt import SignedGraphGT

__all__ = ["SignedGraphGT"]
