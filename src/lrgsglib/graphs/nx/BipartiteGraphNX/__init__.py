"""
BipartiteGraphNX: NetworkX-based signed bipartite (two-mode) graph.

This module provides the BipartiteGraphNX class for generating signed
bipartite random graphs with configurable edge probability.
"""

from .BipartiteGraphNX import BipartiteGraphNX

# Backward compatibility alias
BipartiteGraph = BipartiteGraphNX

__all__ = [
    "BipartiteGraphNX",
    "BipartiteGraph",  # alias for backward compatibility
]
