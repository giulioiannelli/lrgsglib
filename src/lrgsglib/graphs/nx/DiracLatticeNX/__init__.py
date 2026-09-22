"""DiracLatticeNX: NetworkX-based Dirac lattice structures.

This module provides the DiracLatticeGraphNX abstract base class and
concrete implementations DiracCombGraphNX and DiracBrushGraphNX for
hierarchical graphs with efficient spectral computation.
"""

from .DiracBrushNX import DiracBrushGraph, DiracBrushGraphNX
from .DiracCombNX import DiracCombGraph, DiracCombGraphNX
from .DiracLatticeNX import DiracLatticeGraph, DiracLatticeGraphNX

__all__ = [
    # NX classes
    "DiracLatticeGraphNX",
    "DiracCombGraphNX",
    "DiracBrushGraphNX",
    # Backward compatibility aliases
    "DiracLatticeGraph",
    "DiracCombGraph",
    "DiracBrushGraph",
]
