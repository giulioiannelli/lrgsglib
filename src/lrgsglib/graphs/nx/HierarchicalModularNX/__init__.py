"""Hierarchical Modular Network module."""

from .generators_hmn import hierarchical_modular_network
from .HierarchicalModularNetworkNX import HierarchicalModularNetworkNX

# Backward compatibility alias
HierarchicalModularNetwork = HierarchicalModularNetworkNX

__all__ = [
    "HierarchicalModularNetworkNX",
    "HierarchicalModularNetwork",  # alias for backward compatibility
    "hierarchical_modular_network",
]
