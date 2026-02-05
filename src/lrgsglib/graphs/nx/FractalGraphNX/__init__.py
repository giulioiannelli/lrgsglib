"""
FractalGraphNX: Abstract base class for fractal/self-similar signed graphs.

This module provides the FractalGraphNX abstract base class and FractalNwContainerBase
for creating fractal graphs with signed edges.
"""

from .FractalGraphNX import FractalGraphNX, FractalNwContainerBase

# Backward compatibility alias
FractalGraph = FractalGraphNX

__all__ = [
    "FractalGraphNX",
    "FractalGraph",  # alias for backward compatibility
    "FractalNwContainerBase",
]
