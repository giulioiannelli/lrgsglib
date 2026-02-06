"""
DualBarabasiAlbertGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import DualBarabasiAlbertNX

warnings.warn(
    "DualBarabasiAlbertGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

DualBarabasiAlbertGT = DualBarabasiAlbertNX

__all__ = ["DualBarabasiAlbertGT"]
