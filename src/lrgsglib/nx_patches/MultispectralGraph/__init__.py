"""Multispectral graph generators and SignedGraph subclasses.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.MultispectralGraph is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

# Re-export from graphs.nx to avoid circular imports
from ...graphs.nx.MultispectralGraphNX import MultispectralGraphNX as MultispectralGraph
from ...graphs.nx.MultiplicativeCascadeNX import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
from ...graphs.nx.HierarchicalModularNX import HierarchicalModularNetworkNX as HierarchicalModularNetwork
from ...graphs.nx.VicsekNX import VicsekGraphNX as VicsekGraph
from ...graphs.nx.DiracLatticeNX import (
    DiracLatticeGraphNX as DiracLatticeGraph,
    DiracCombGraphNX as DiracCombGraph,
    DiracBrushGraphNX as DiracBrushGraph,
)

# Also export generator functions for backward compatibility
from ...graphs.nx.MultispectralGraphNX.generators_msg import (
    multiplicative_cascade_probability_matrix,
    multiplicative_cascade_graph,
    dirac_comb_graph,
    dirac_brush_graph,
)
from ...graphs.nx.HierarchicalModularNX.generators_hmn import hierarchical_modular_network

__all__ = [
    # Main classes
    "MultispectralGraph",
    "MultiplicativeCascadeGraph",
    "HierarchicalModularNetwork",
    "VicsekGraph",
    "DiracLatticeGraph",
    "DiracCombGraph",
    "DiracBrushGraph",
    # Generator functions
    "multiplicative_cascade_probability_matrix",
    "multiplicative_cascade_graph",
    "hierarchical_modular_network",
    "dirac_comb_graph",
    "dirac_brush_graph",
]
