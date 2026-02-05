"""
ExtendedBarabasiAlbertGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.random import ExtendedBarabasiAlbertNX

warnings.warn(
    "ExtendedBarabasiAlbertGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

ExtendedBarabasiAlbertGT = ExtendedBarabasiAlbertNX

__all__ = ["ExtendedBarabasiAlbertGT"]
