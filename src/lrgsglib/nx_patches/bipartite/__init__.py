"""
Bipartite graph family module.

This module provides signed bipartite (two-mode) networks where nodes
are divided into two disjoint sets and edges only connect nodes from
different sets.

Graph Types
-----------
- BipartiteGraph: Random bipartite graph with edge probability p
- BipartiteFromDegreeSequence: Bipartite with prescribed degree sequences

Example
-------
>>> from lrgsglib.nx_patches.bipartite import BipartiteGraph
>>> bg = BipartiteGraph(n1=50, n2=100, p=0.05, pflip=0.1, seed=42)
>>> bg.flip_random_fract_edges()
>>> B = bg.get_biadjacency_matrix()  # Shape: (50, 100)
"""

from .bipartite_graph import BipartiteGraph, BipartiteFromDegreeSequence

__all__ = [
    "BipartiteGraph",
    "BipartiteFromDegreeSequence",
]
