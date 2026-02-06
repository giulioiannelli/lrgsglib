"""
ErdosRenyiGT: graph-tool-based Erdos-Renyi signed graph.

This module re-exports the ErdosRenyiGT class from gt_patches for the
unified graph module structure.

The actual implementation remains in gt_patches.random.erdos_renyi_gt for backward
compatibility.

Example
-------
>>> from lrgsglib.graphs.gt.random import ErdosRenyiGT
>>> er = ErdosRenyiGT(n=200, p=0.05, pflip=0.2, seed=42)
>>> er.flip_random_fract_edges()
>>> print(f"Nodes: {er.N}, Edges: {er.num_edges}")
"""

from ....gt_patches.random.erdos_renyi_gt import ErdosRenyiGT

__all__ = ["ErdosRenyiGT"]
