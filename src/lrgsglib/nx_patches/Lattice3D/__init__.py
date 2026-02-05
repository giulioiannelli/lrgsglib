"""Lattice3D module.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.Lattice3D is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

from ...graphs.nx.Lattice3DNX import (
    Lattice3DNX as Lattice3D,
    load_or_compute_Lattice3DNX as load_or_compute_Lattice3D,
)

__all__ = [
    "Lattice3D",
    "load_or_compute_Lattice3D",
]
