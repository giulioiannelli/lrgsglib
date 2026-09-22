"""
Lattice2DNX: NetworkX-based 2D lattice signed graph.

This module provides the Lattice2DNX class and supporting functions for
2D lattice graphs with periodic/fixed boundary conditions.
"""

from . import generators_2d, paths
from .eigenspace import (
    create_lattice_with_eigenspace,
    load_or_compute_Lattice2DNX,
)
from .Lattice2DNX import Lattice2DNX

# Backward compatibility alias
Lattice2D = Lattice2DNX
load_or_compute_Lattice2D = load_or_compute_Lattice2DNX

__all__ = [
    "Lattice2DNX",
    "Lattice2D",  # alias for backward compatibility
    "create_lattice_with_eigenspace",
    "load_or_compute_Lattice2DNX",
    "load_or_compute_Lattice2D",  # alias
    "generators_2d",
    "paths",
]
