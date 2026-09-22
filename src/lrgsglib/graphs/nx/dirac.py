"""Dirac lattice graph facade module.

Re-exports Dirac lattice classes from DiracLatticeNX package for cleaner imports.
"""

from .DiracLatticeNX import (  # Backward compatibility aliases
    DiracBrushGraph,
    DiracBrushGraphNX,
    DiracCombGraph,
    DiracCombGraphNX,
    DiracLatticeGraph,
    DiracLatticeGraphNX,
)

__all__ = [
    "DiracLatticeGraphNX",
    "DiracCombGraphNX",
    "DiracBrushGraphNX",
    "DiracLatticeGraph",
    "DiracCombGraph",
    "DiracBrushGraph",
]
