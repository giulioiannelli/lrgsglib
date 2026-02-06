"""
ExtendedBarabasiAlbertNX: NetworkX-based extended BA model.

This module provides the ExtendedBarabasiAlbertNX class for generating signed
scale-free networks with initial attractiveness parameter.
"""

from .ExtendedBarabasiAlbertNX import ExtendedBarabasiAlbertNX

# Backward compatibility alias
ExtendedBarabasiAlbert = ExtendedBarabasiAlbertNX

__all__ = [
    "ExtendedBarabasiAlbertNX",
    "ExtendedBarabasiAlbert",  # alias for backward compatibility
]
