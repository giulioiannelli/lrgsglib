"""Lattice2D module.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import sys
import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.Lattice2D is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

from ...graphs.nx.Lattice2DNX import (
    Lattice2DNX as Lattice2D,
    load_or_compute_Lattice2DNX as load_or_compute_Lattice2D,
    create_lattice_with_eigenspace,
)

__all__ = [
    "Lattice2D",
    "load_or_compute_Lattice2D",
    "create_lattice_with_eigenspace",
]

# Fix submodule shadowing
_parent = sys.modules.get("lrgsglib.nx_patches")
if _parent is not None:
    _parent.Lattice2D = Lattice2D
