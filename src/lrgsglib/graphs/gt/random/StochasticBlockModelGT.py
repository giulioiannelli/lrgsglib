"""
StochasticBlockModelGT: graph-tool-based Stochastic Block Model signed graph.

This module re-exports the StochasticBlockModelGT class from gt_patches for the
unified graph module structure.

The actual implementation remains in gt_patches.random.stochastic_block_gt for backward
compatibility.

Example
-------
>>> from lrgsglib.graphs.gt.random import StochasticBlockModelGT
>>> sizes = [30, 30, 40]
>>> p_matrix = [[0.3, 0.05, 0.05], [0.05, 0.3, 0.05], [0.05, 0.05, 0.3]]
>>> sbm = StochasticBlockModelGT(sizes=sizes, p_matrix=p_matrix, seed=42)
>>> sbm.flip_random_fract_edges()
>>> print(f"Nodes: {sbm.N}, Communities: {sbm.num_communities}")
"""

from ....gt_patches.random.stochastic_block_gt import StochasticBlockModelGT

__all__ = ["StochasticBlockModelGT"]
