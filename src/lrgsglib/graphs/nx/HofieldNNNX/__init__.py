"""
HofieldNNNX: NetworkX-based Hopfield neural network signed graph.

This module provides the HofieldNNNX class for creating Hopfield-like
neural networks built on top of fully connected graphs.
"""

from .HofieldNNNX import HofieldNNNX
from .patterns import init_mnist_patterns

# Backward compatibility alias
HofieldNN = HofieldNNNX

__all__ = [
    "HofieldNNNX",
    "HofieldNN",  # alias for backward compatibility
    "init_mnist_patterns",
]
