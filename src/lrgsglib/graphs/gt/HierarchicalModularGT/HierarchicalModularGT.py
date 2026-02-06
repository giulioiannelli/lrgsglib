"""
HierarchicalModularNetworkGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.multispectral import HierarchicalModularNetworkNX

warnings.warn(
    "HierarchicalModularNetworkGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

HierarchicalModularNetworkGT = HierarchicalModularNetworkNX

__all__ = ["HierarchicalModularNetworkGT"]
