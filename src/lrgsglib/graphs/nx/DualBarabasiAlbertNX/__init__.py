"""
DualBarabasiAlbertNX: NetworkX-based dual BA model.

This module provides the DualBarabasiAlbertNX class for generating signed
scale-free networks with two attachment modes.
"""

from .DualBarabasiAlbertNX import DualBarabasiAlbertNX

# Backward compatibility alias
DualBarabasiAlbert = DualBarabasiAlbertNX

__all__ = [
    "DualBarabasiAlbertNX",
    "DualBarabasiAlbert",  # alias for backward compatibility
]
