"""
CompleteGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.complete import CompleteGraphNX

warnings.warn(
    "CompleteGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

CompleteGraphGT = CompleteGraphNX

__all__ = ["CompleteGraphGT"]
