"""Multispectral graph generators and SignedGraph subclasses."""

from .MultispectralGraph import MultispectralGraph
from ..MultiplicativeCascade import MultiplicativeCascadeGraph
from ..HierarchicalModular import HierarchicalModularNetwork
from ..Vicsek import VicsekGraph
from ..DiracLattice import DiracLatticeGraph, DiracCombGraph, DiracBrushGraph

# Also export generator functions for backward compatibility
from .generators_msg import (
    multiplicative_cascade_probability_matrix,
    multiplicative_cascade_graph,
    dirac_comb_graph,
    dirac_brush_graph,
)
from ..HierarchicalModular.generators_hmn import hierarchical_modular_network

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
