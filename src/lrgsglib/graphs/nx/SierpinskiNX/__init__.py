"""
SierpinskiNX: NetworkX-based Sierpinski gasket/triangle fractal signed graph.

This module provides the SierpinskiNX class for creating Sierpinski fractal
graphs (gasket, carpet, tetrahedron variants) with signed edges.
"""

from .SierpinskiNX import SierpinskiNX

# Backward compatibility aliases
SierpinskiGraph = SierpinskiNX
SierpinskiGraphNX = SierpinskiNX  # Alternate NX naming pattern

__all__ = [
    "SierpinskiNX",
    "SierpinskiGraph",  # alias for backward compatibility
    "SierpinskiGraphNX",  # alternate naming pattern
]
