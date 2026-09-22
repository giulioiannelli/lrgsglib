"""
FullyConnectedNX: NetworkX-based fully connected (complete) signed graph.

This module provides the FullyConnectedNX class for creating fully connected
graphs with signed edges and animation capabilities.
"""

from .animation import make_animation
from .FullyConnectedNX import FullyConnectedNX

# Backward compatibility alias
FullyConnected = FullyConnectedNX

__all__ = [
    "FullyConnectedNX",
    "FullyConnected",  # alias for backward compatibility
    "make_animation",
]
