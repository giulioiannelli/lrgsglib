"""
BarabasiAlbertNX: NetworkX-based Barabasi-Albert scale-free graph.

This module provides the BarabasiAlbertNX class for generating signed
scale-free networks using preferential attachment.
"""

from .BarabasiAlbertNX import BarabasiAlbertNX

# Backward compatibility alias
BarabasiAlbert = BarabasiAlbertNX

__all__ = [
    "BarabasiAlbertNX",
    "BarabasiAlbert",  # alias for backward compatibility
]
