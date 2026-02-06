"""
MultispectralGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.multispectral import MultispectralGraphNX

warnings.warn(
    "MultispectralGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

MultispectralGraphGT = MultispectralGraphNX

__all__ = ["MultispectralGraphGT"]
