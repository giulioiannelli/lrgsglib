"""
LFRBenchmarkGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import LFRBenchmarkNX

warnings.warn(
    "LFRBenchmarkGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

LFRBenchmarkGT = LFRBenchmarkNX

__all__ = ["LFRBenchmarkGT"]
