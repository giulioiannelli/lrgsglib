"""
graph-tool C++ extensions for lrgsglib.

This subpackage provides high-performance C++ extensions for graph-tool operations.
Extensions are built via `make cpp-make` from the lrgsglib root directory.

Available functions:
- create_triangular_lattice(width, height) - Create a triangular lattice graph

The C++ extensions require graph-tool to be installed and provide significant
speedups (typically 10-50x) over pure Python implementations.
"""
from __future__ import annotations

# Import available extensions
try:
    from .lattice_generators import create_triangular_lattice
except ImportError as e:
    # Extension not built - provide helpful error
    def create_triangular_lattice(*args, **kwargs):
        raise ImportError(
            "triangular_lattice C++ extension not available. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from e


__all__ = [
    "create_triangular_lattice",
]
