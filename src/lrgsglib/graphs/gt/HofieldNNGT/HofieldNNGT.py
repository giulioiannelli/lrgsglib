"""
HofieldNNGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.neural import HofieldNNNX

warnings.warn(
    "HofieldNNGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

HofieldNNGT = HofieldNNNX

__all__ = ["HofieldNNGT"]
