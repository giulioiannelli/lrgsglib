"""
Triangular lattice generator using graph-tool C++ extension.

This module provides a Python wrapper for the C++ triangular lattice generator.
The C++ implementation is significantly faster than pure Python alternatives.
"""
from __future__ import annotations

import sys
from pathlib import Path

import graph_tool.all as gt


def _load_extension():
    """Load the triangular_lattice C++ extension module."""
    # The .so file is in the same directory as this Python file
    cpp_dir = Path(__file__).parent
    so_file = cpp_dir / "triangular_lattice.so"

    if not so_file.exists():
        raise ImportError(
            f"triangular_lattice.so not found at {so_file}. "
            f"Run 'make cpp-make' from the lrgsglib root directory to build it."
        )

    # Add the cpp directory to sys.path if not already there
    cpp_dir_str = str(cpp_dir)
    if cpp_dir_str not in sys.path:
        sys.path.insert(0, cpp_dir_str)

    # Now import the module
    import triangular_lattice as _tl_module

    return _tl_module


# Load the extension at module import time
_triangular_lattice_ext = _load_extension()


def create_triangular_lattice(width: int, height: int) -> gt.Graph:
    """
    Create a triangular lattice graph with given width and height.

    Uses a C++ extension for efficient implementation (~50x faster than pure Python).

    Args:
        width: The width of the lattice (number of columns).
        height: The height of the lattice (number of rows).

    Returns:
        A graph_tool Graph object representing the triangular lattice.
        The graph is undirected with width*height vertices.

    Example:
        >>> g = create_triangular_lattice(10, 10)
        >>> g.num_vertices()
        100
    """
    # Create a new empty Graph using graph_tool
    g = gt.Graph(directed=False)

    # For graph objects, we need to pass the internal GraphInterface
    # which is accessed via g._Graph__graph.
    # The C++ function needs no property maps, so we only pass dimensions
    _triangular_lattice_ext.create_triangular_lattice(g._Graph__graph, width, height)

    return g
