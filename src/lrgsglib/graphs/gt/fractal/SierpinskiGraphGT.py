"""
SierpinskiGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.fractal import SierpinskiGraphNX

warnings.warn(
    "SierpinskiGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

SierpinskiGraphGT = SierpinskiGraphNX

__all__ = ["SierpinskiGraphGT"]
