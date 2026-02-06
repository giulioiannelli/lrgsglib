"""
RandomGeometricNX: NetworkX-based random geometric graph.

This module provides the RandomGeometricNX class for generating signed
random graphs with spatial embedding.
"""

from .RandomGeometricNX import RandomGeometricNX

# Backward compatibility alias
RandomGeometric = RandomGeometricNX

__all__ = [
    "RandomGeometricNX",
    "RandomGeometric",  # alias for backward compatibility
]
