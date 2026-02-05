"""
RandomGraphNX: NetworkX-based abstract random graph base class.

This module provides the RandomGraphNX abstract base class and
RandomNwContainerBase helper class for random graph models.
"""

from .RandomGraphNX import RandomGraphNX, RandomNwContainerBase

# Backward compatibility alias
RandomGraph = RandomGraphNX

__all__ = [
    "RandomGraphNX",
    "RandomGraph",  # alias for backward compatibility
    "RandomNwContainerBase",
]
