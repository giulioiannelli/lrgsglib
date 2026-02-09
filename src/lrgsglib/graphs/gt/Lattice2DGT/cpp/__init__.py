"""
C++ extensions for graph-tool lattice generation.

Available functions:
- create_triangular_lattice(width, height) - Create a triangular lattice graph
"""
from __future__ import annotations

try:
    from .lattice_generators import create_triangular_lattice
except ImportError as e:
    def create_triangular_lattice(*args, **kwargs):
        raise ImportError(
            "triangular_lattice C++ extension not available. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from e


__all__ = [
    "create_triangular_lattice",
]
