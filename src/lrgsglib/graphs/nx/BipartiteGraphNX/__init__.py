"""
BipartiteGraphNX: NetworkX-based signed bipartite (two-mode) graph.

This module provides the BipartiteGraphNX class for generating signed
bipartite graphs with configurable edge probability or prescribed
degree sequences.
"""

from .BipartiteGraphNX import BipartiteGraphNX

# Backward compatibility aliases
BipartiteGraph = BipartiteGraphNX
BipartiteFromDegreeSequenceNX = BipartiteGraphNX
BipartiteFromDegreeSequence = BipartiteGraphNX

__all__ = [
    "BipartiteGraphNX",
    "BipartiteGraph",
    "BipartiteFromDegreeSequenceNX",
    "BipartiteFromDegreeSequence",
]
