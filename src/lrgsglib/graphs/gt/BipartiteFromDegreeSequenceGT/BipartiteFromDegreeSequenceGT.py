"""
BipartiteFromDegreeSequenceGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.bipartite import BipartiteFromDegreeSequenceNX

warnings.warn(
    "BipartiteFromDegreeSequenceGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

BipartiteFromDegreeSequenceGT = BipartiteFromDegreeSequenceNX

__all__ = ["BipartiteFromDegreeSequenceGT"]
