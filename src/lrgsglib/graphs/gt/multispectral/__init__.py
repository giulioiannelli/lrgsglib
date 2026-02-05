"""
Native graph-tool multispectral/hierarchical graph implementations.

Includes:
- VicsekGraphGT: Native GT implementation
- MultiplicativeCascadeGraphGT: Native GT implementation
- HierarchicalModularNetworkGT: Falls back to NetworkX
- MultispectralGraphGT: Abstract base class
"""

from .MultispectralGraphGT import MultispectralGraphGT
from .MultiplicativeCascadeGraphGT import MultiplicativeCascadeGraphGT
from .VicsekGraphGT import VicsekGraphGT
from .HierarchicalModularNetworkGT import HierarchicalModularNetworkGT

__all__ = [
    "MultispectralGraphGT",
    "MultiplicativeCascadeGraphGT",
    "VicsekGraphGT",
    "HierarchicalModularNetworkGT",
]
