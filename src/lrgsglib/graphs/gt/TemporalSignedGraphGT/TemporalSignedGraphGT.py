"""
TemporalSignedGraphGT: Deprecated - use TemporalGraphGT directly.

TemporalGraphNX now inherits from SignedGraphNX, so the separate
TemporalSignedGraph class is no longer needed.
"""

import warnings
from ...nx.TemporalGraphNX import TemporalGraphNX

warnings.warn(
    "TemporalSignedGraphGT is deprecated. Use TemporalGraphGT (or TemporalGraphNX) "
    "directly — it now inherits from SignedGraphNX.",
    DeprecationWarning,
    stacklevel=2,
)

TemporalSignedGraphGT = TemporalGraphNX

__all__ = ["TemporalSignedGraphGT"]
