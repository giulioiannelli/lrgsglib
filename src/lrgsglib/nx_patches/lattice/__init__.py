"""
Unified lattice graph module.

This module provides a clean class hierarchy for lattice graphs:
- LatticeND: Abstract base class for N-dimensional lattices
- Lattice2D: Backward-compatible 2D lattice wrapper
- Lattice3D: Backward-compatible 3D lattice wrapper

Example
-------
>>> from lrgsglib.nx_patches.lattice import Lattice2D, Lattice3D
>>> lat2d = Lattice2D(side1=10, geo='sqr', pflip=0.1)
>>> lat3d = Lattice3D(dim=8, geo='sc', pflip=0.1)
"""

from .lattice_nd import LatticeND
from .lattice_2d import Lattice2D
from .lattice_3d import Lattice3D

__all__ = [
    "LatticeND",
    "Lattice2D",
    "Lattice3D",
]
