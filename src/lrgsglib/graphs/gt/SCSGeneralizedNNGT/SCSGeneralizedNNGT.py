"""
SCSGeneralizedNNGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.neural import SCSGeneralizedNNNX

warnings.warn(
    "SCSGeneralizedNNGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

SCSGeneralizedNNGT = SCSGeneralizedNNNX

__all__ = ["SCSGeneralizedNNGT"]
