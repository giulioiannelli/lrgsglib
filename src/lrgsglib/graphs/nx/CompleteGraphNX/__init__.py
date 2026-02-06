"""
CompleteGraphNX: Abstract base class for fully connected signed graphs.

This module provides the CompleteGraphNX abstract base class for creating
complete (fully connected) graphs with signed edges.
"""

from .CompleteGraphNX import CompleteGraphNX

# Backward compatibility alias
CompleteGraph = CompleteGraphNX

__all__ = [
    "CompleteGraphNX",
    "CompleteGraph",  # alias for backward compatibility
]
