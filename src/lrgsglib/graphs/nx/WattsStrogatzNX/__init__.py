"""
WattsStrogatzNX: NetworkX-based Watts-Strogatz small-world signed graph.

This module provides the WattsStrogatzNX class for generating signed
Watts-Strogatz small-world graphs.
"""

from .WattsStrogatzNX import WattsStrogatzNX

# Backward compatibility alias
WattsStrogatz = WattsStrogatzNX

__all__ = [
    "WattsStrogatzNX",
    "WattsStrogatz",  # alias for backward compatibility
]
