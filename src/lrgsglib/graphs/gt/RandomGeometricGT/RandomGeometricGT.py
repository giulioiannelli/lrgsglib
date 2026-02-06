"""
RandomGeometricGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import RandomGeometricNX

warnings.warn(
    "RandomGeometricGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

RandomGeometricGT = RandomGeometricNX

__all__ = ["RandomGeometricGT"]
