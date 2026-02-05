"""
C++ extensions for graph-tool lattice generation.

This module re-exports C++ extensions from gt_patches.cpp for the
unified graph module structure.

The actual implementations remain in gt_patches.cpp for backward compatibility.

Available functions:
- create_triangular_lattice(width, height) - Create a triangular lattice graph
- kcore_decomposition(g, core=None) - K-core decomposition
"""
from __future__ import annotations

# Re-export from the original location
try:
    from .....gt_patches.cpp.lattice_generators import create_triangular_lattice
except ImportError as e:
    def create_triangular_lattice(*args, **kwargs):
        raise ImportError(
            "triangular_lattice C++ extension not available. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from e

try:
    from .....gt_patches.cpp.kcore import kcore_decomposition
except ImportError as e:
    def kcore_decomposition(*args, **kwargs):
        raise ImportError(
            "kcore C++ extension not available. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from e


__all__ = [
    "create_triangular_lattice",
    "kcore_decomposition",
]
