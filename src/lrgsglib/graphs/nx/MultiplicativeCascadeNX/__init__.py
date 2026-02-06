"""MultiplicativeCascadeNX: NetworkX-based multiplicative cascade graph.

This module provides the MultiplicativeCascadeGraphNX class for generating
signed multiplicative cascade graphs with hierarchical structure.
"""

from .MultiplicativeCascadeNX import MultiplicativeCascadeGraphNX

# Backward compatibility alias
MultiplicativeCascadeGraph = MultiplicativeCascadeGraphNX

__all__ = [
    "MultiplicativeCascadeGraphNX",
    "MultiplicativeCascadeGraph",  # alias for backward compatibility
]
