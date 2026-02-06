"""
HolmeKimGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import HolmeKimNX

warnings.warn(
    "HolmeKimGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

HolmeKimGT = HolmeKimNX

__all__ = ["HolmeKimGT"]
