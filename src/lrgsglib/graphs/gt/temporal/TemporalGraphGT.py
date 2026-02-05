"""
TemporalGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.temporal import TemporalGraphNX

warnings.warn(
    "TemporalGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

TemporalGraphGT = TemporalGraphNX

__all__ = ["TemporalGraphGT"]
