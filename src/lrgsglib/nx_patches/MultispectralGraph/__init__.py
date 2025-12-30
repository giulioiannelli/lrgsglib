"""Multispectral graph generators and SignedGraph subclasses."""

from .MultispectralGraph import MultispectralGraph
from ..MultiplicativeCascade import MultiplicativeCascadeGraph
from ..Vicsek import VicsekGraph
from ..DiracLattice import DiracLatticeGraph, DiracCombGraph, DiracBrushGraph

# Also export generator functions for backward compatibility
from .generators_msg import (
    multiplicative_cascade_probability_matrix,
    multiplicative_cascade_graph,
    dirac_comb_graph,
    dirac_brush_graph,
)

__all__ = [
    # Main classes
    "MultispectralGraph",
    "MultiplicativeCascadeGraph",
    "VicsekGraph",
    "DiracLatticeGraph",
    "DiracCombGraph",
    "DiracBrushGraph",
    # Generator functions
    "multiplicative_cascade_probability_matrix",
    "multiplicative_cascade_graph",
    "dirac_comb_graph",
    "dirac_brush_graph",
]
