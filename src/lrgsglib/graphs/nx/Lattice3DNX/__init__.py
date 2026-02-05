"""
Lattice3DNX: NetworkX-based 3D lattice signed graph.

This module provides the Lattice3DNX class and supporting functions for
3D lattice graphs (SC, BCC, FCC) with periodic/fixed boundary conditions.
"""

from .Lattice3DNX import Lattice3DNX
from .eigenspace import load_or_compute_Lattice3DNX
from . import generators_3d

# Backward compatibility alias
Lattice3D = Lattice3DNX
load_or_compute_Lattice3D = load_or_compute_Lattice3DNX

__all__ = [
    "Lattice3DNX",
    "Lattice3D",  # alias
    "load_or_compute_Lattice3DNX",
    "load_or_compute_Lattice3D",  # alias
    "generators_3d",
]
