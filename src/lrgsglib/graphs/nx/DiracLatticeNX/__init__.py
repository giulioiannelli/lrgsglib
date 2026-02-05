"""DiracLatticeNX: NetworkX-based Dirac lattice structures.

This module provides the DiracLatticeGraphNX abstract base class and
concrete implementations DiracCombGraphNX and DiracBrushGraphNX for
hierarchical graphs with efficient spectral computation.
"""

from .DiracLatticeNX import DiracLatticeGraphNX, DiracLatticeGraph
from .DiracCombNX import DiracCombGraphNX, DiracCombGraph
from .DiracBrushNX import DiracBrushGraphNX, DiracBrushGraph

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
