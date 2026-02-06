"""
FullyConnectedNX: NetworkX-based fully connected (complete) signed graph.

This module provides the FullyConnectedNX class for creating fully connected
graphs with signed edges and animation capabilities.
"""

from .FullyConnectedNX import FullyConnectedNX
from .animation import make_animation

# Backward compatibility alias
FullyConnected = FullyConnectedNX

__all__ = [
    "FullyConnectedNX",
    "FullyConnected",  # alias for backward compatibility
    "make_animation",
]
