"""Lattice graph facade module.

Re-exports lattice graph classes from Lattice2DNX and Lattice3DNX packages
for cleaner imports.
"""

from .Lattice2DNX import (
    Lattice2DNX,
    # Backward compatibility alias
    Lattice2D,
    create_lattice_with_eigenspace,
    load_or_compute_Lattice2DNX,
    load_or_compute_Lattice2D,
)
from .Lattice3DNX import (
    Lattice3DNX,
    # Backward compatibility alias
    Lattice3D,
    load_or_compute_Lattice3DNX,
    load_or_compute_Lattice3D,
)

# LatticeNDNX is referenced but not implemented yet - create a placeholder
# Until the abstract base class is implemented, point to Lattice2DNX as a workaround
LatticeNDNX = Lattice2DNX

__all__ = [
    # Abstract base (placeholder)
    "LatticeNDNX",
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
