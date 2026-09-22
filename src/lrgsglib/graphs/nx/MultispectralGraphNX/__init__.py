"""MultispectralGraphNX: Abstract base class for multispectral graph generators.

This module provides the MultispectralGraphNX class and supporting generator functions
for graphs with multiple spectral dimensions.
"""

from .generators_msg import (
    dirac_brush_graph,
    dirac_comb_graph,
    initial_measure,
    link_probabilities,
    multiplicative_cascade_exp_clocks,
    multiplicative_cascade_graph,
    multiplicative_cascade_probability_matrix,
    palla_lovasz_vicksek_graph,
)
from .MultispectralGraphNX import MultispectralGraphNX

# Backward compatibility alias
MultispectralGraph = MultispectralGraphNX

__all__ = [
    # Main class
    "MultispectralGraphNX",
    "MultispectralGraph",  # alias
    # Generator functions
    "multiplicative_cascade_probability_matrix",
    "multiplicative_cascade_graph",
    "multiplicative_cascade_exp_clocks",
    "initial_measure",
    "link_probabilities",
    "palla_lovasz_vicksek_graph",
    "dirac_comb_graph",
    "dirac_brush_graph",
]
