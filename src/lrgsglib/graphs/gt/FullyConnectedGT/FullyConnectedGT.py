"""
FullyConnectedGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.complete import FullyConnectedNX

warnings.warn(
    "FullyConnectedGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

FullyConnectedGT = FullyConnectedNX

__all__ = ["FullyConnectedGT"]
