"""Dirac lattice graph facade module.

Re-exports Dirac lattice classes from DiracLatticeNX package for cleaner imports.
"""

from .DiracLatticeNX import (
    DiracLatticeGraphNX,
    DiracCombGraphNX,
    DiracBrushGraphNX,
    # Backward compatibility aliases
    DiracLatticeGraph,
    DiracCombGraph,
    DiracBrushGraph,
)

__all__ = [
    "DiracLatticeGraphNX",
    "DiracCombGraphNX",
    "DiracBrushGraphNX",
    "DiracLatticeGraph",
    "DiracCombGraph",
    "DiracBrushGraph",
]
