"""
kRegularGraphNX: NetworkX-based k-regular random graph.

This module provides the kRegularGraphNX class for generating signed
random graphs where every node has exactly degree k.
"""

from .kRegularGraphNX import kRegularGraphNX

# Backward compatibility alias
kRegularGraph = kRegularGraphNX

__all__ = [
    "kRegularGraphNX",
    "kRegularGraph",  # alias for backward compatibility
]
