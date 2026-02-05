"""
graph-tool bipartite graph implementations with GT suffix.

Note: These fall back to NetworkX implementations.
"""

from .BipartiteGraphGT import BipartiteGraphGT
from .BipartiteFromDegreeSequenceGT import BipartiteFromDegreeSequenceGT

__all__ = [
    "BipartiteGraphGT",
    "BipartiteFromDegreeSequenceGT",
]
