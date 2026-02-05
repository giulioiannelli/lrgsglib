"""Hierarchical Modular Network module.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.HierarchicalModular is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

from ...graphs.nx.HierarchicalModularNX import HierarchicalModularNetworkNX as HierarchicalModularNetwork
from ...graphs.nx.HierarchicalModularNX.generators_hmn import hierarchical_modular_network

__all__ = [
    "HierarchicalModularNetwork",
    "hierarchical_modular_network",
]
