"""
TemporalSignedGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.temporal import TemporalSignedGraphNX

warnings.warn(
    "TemporalSignedGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

TemporalSignedGraphGT = TemporalSignedGraphNX

__all__ = ["TemporalSignedGraphGT"]
