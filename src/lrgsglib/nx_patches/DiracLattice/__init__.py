"""Dirac lattice graph module.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.DiracLattice is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

from ...graphs.nx.DiracLatticeNX import (
    DiracLatticeGraphNX as DiracLatticeGraph,
    DiracCombGraphNX as DiracCombGraph,
    DiracBrushGraphNX as DiracBrushGraph,
)

__all__ = ["DiracLatticeGraph", "DiracCombGraph", "DiracBrushGraph"]
