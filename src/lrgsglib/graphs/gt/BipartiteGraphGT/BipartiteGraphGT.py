"""
BipartiteGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.bipartite import BipartiteGraphNX

warnings.warn(
    "BipartiteGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

BipartiteGraphGT = BipartiteGraphNX

__all__ = ["BipartiteGraphGT"]
