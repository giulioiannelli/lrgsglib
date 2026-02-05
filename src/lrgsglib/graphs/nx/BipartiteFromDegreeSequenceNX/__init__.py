"""
BipartiteFromDegreeSequenceNX: NetworkX-based bipartite graph with prescribed degrees.

This module provides the BipartiteFromDegreeSequenceNX class for generating
signed bipartite graphs with specified degree sequences.
"""

from .BipartiteFromDegreeSequenceNX import BipartiteFromDegreeSequenceNX

# Backward compatibility alias
BipartiteFromDegreeSequence = BipartiteFromDegreeSequenceNX

__all__ = [
    "BipartiteFromDegreeSequenceNX",
    "BipartiteFromDegreeSequence",  # alias for backward compatibility
]
