"""
DGMgraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.fractal import DGMgraphNX

warnings.warn(
    "DGMgraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

DGMgraphGT = DGMgraphNX

__all__ = ["DGMgraphGT"]
