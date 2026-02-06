"""Hierarchical Modular Network module."""

from .HierarchicalModularNetworkNX import HierarchicalModularNetworkNX
from .generators_hmn import hierarchical_modular_network

# Backward compatibility alias
HierarchicalModularNetwork = HierarchicalModularNetworkNX

__all__ = [
    "HierarchicalModularNetworkNX",
    "HierarchicalModularNetwork",  # alias for backward compatibility
    "hierarchical_modular_network",
]
