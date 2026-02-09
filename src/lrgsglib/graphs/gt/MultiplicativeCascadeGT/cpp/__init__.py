"""
C++ extension for MultiplicativeCascadeGT edge building.

Provides high-performance grid-edge construction using graph-tool's
C++ backend.  Falls back to vectorized numpy if the .so is not built.

Build with: make cpp-make (from lrgsglib root)
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import graph_tool

_SO_PATH = Path(__file__).parent / "multiplicative_cascade.so"

try:
    from . import multiplicative_cascade as _mc_module
    _build_cascade_edges_raw = _mc_module.build_cascade_edges
except ImportError as e:
    _mc_module = None
    _import_error = e

    def _build_cascade_edges_raw(*args, **kwargs):
        raise ImportError(
            f"multiplicative_cascade C++ extension not available at "
            f"{_SO_PATH}. Build with 'make cpp-make' from the lrgsglib "
            "root directory."
        ) from _import_error


def build_cascade_edges(
    g: "graph_tool.Graph",
    rows: list,
    cols: list,
    n_selected: int,
    size: int,
    periodic: bool = False,
) -> None:
    """Build grid edges between selected nodes on a 2D lattice.

    Parameters
    ----------
    g : graph_tool.Graph
        Graph with ``n_selected`` vertices already added.
    rows : list[int]
        Row coordinates of selected nodes.
    cols : list[int]
        Column coordinates of selected nodes.
    n_selected : int
        Number of selected nodes.
    size : int
        Grid side length.
    periodic : bool
        Use periodic boundary conditions.
    """
    _build_cascade_edges_raw(
        g._Graph__graph, rows, cols, n_selected, size, periodic
    )


__all__ = ["build_cascade_edges"]
