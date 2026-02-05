"""
RandomGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import RandomGraphNX

warnings.warn(
    "RandomGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

RandomGraphGT = RandomGraphNX

__all__ = ["RandomGraphGT"]
