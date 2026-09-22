"""Lattice graph facade module.

Re-exports lattice graph classes from Lattice2DNX and Lattice3DNX packages
for cleaner imports.
"""

from .Lattice2DNX import (  # Backward compatibility alias
    Lattice2D,
    Lattice2DNX,
    create_lattice_with_eigenspace,
    load_or_compute_Lattice2D,
    load_or_compute_Lattice2DNX,
)
from .Lattice3DNX import (  # Backward compatibility alias
    Lattice3D,
    Lattice3DNX,
    load_or_compute_Lattice3D,
    load_or_compute_Lattice3DNX,
)

__all__ = [
    # 2D
    "Lattice2DNX",
    "Lattice2D",
    "create_lattice_with_eigenspace",
    "load_or_compute_Lattice2DNX",
    "load_or_compute_Lattice2D",
    # 3D
    "Lattice3DNX",
    "Lattice3D",
    "load_or_compute_Lattice3DNX",
    "load_or_compute_Lattice3D",
]
