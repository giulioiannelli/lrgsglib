"""
kRegularGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import kRegularGraphNX

warnings.warn(
    "kRegularGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

kRegularGraphGT = kRegularGraphNX

__all__ = ["kRegularGraphGT"]
