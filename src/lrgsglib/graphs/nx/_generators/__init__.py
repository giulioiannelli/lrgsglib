"""
Consolidated graph generators for all graph families.

This module provides a unified interface to all graph generation functions
used across the nx_patches module, organized by graph family.

Generators are organized by family:
- lattice: 2D, 3D, and ND lattice generators
- random: Random graph models (ER, WS, BA, etc.)
- hierarchical: Multispectral, Vicsek, HMN, Dirac
- fractal: DGM and other fractal structures

Example
-------
>>> from lrgsglib.nx_patches._generators import lattice
>>> G = lattice.squared_lattice_graph_FastPatch(10, 10, periodic=True)
"""

from . import lattice
from . import random

__all__ = [
    "lattice",
    "random",
]
