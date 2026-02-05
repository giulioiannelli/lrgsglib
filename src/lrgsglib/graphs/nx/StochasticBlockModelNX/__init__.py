"""
StochasticBlockModelNX: NetworkX-based Stochastic Block Model graph.

This module provides the StochasticBlockModelNX class for generating signed
random graphs with planted community structure.
"""

from .StochasticBlockModelNX import StochasticBlockModelNX

# Backward compatibility alias
StochasticBlockModel = StochasticBlockModelNX

__all__ = [
    "StochasticBlockModelNX",
    "StochasticBlockModel",  # alias for backward compatibility
]
