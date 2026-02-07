"""
C++ extension for BarabasiAlbertGT graph generation.

Provides high-performance Barabasi-Albert preferential attachment
graph construction using graph-tool's C++ backend.

Build with: make cpp-make (from lrgsglib root)
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import graph_tool

# Import the compiled extension
_SO_PATH = Path(__file__).parent / "barabasi_albert.so"

try:
    from . import barabasi_albert as _ba_module
    _create_barabasi_albert_raw = _ba_module.create_barabasi_albert
except ImportError as e:
    _ba_module = None
    _import_error = e

    def _create_barabasi_albert_raw(*args, **kwargs):
        raise ImportError(
            f"barabasi_albert C++ extension not available at {_SO_PATH}. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from _import_error


def create_barabasi_albert(
    g: "graph_tool.Graph",
    n: int,
    m: int,
    seed: int = 0,
) -> None:
    """
    Create a Barabasi-Albert preferential attachment graph.

    This is a high-performance C++ implementation providing ~100x speedup
    over the pure Python graph-tool implementation.

    Parameters
    ----------
    g : graph_tool.Graph
        An empty graph-tool Graph object to populate.
    n : int
        Total number of nodes.
    m : int
        Number of edges each new node creates (must be < n).
    seed : int, optional
        Random seed for reproducibility. 0 = use random device.

    Examples
    --------
    >>> import graph_tool.all as gt
    >>> from lrgsglib.graphs.gt.BarabasiAlbertGT.cpp import create_barabasi_albert
    >>> g = gt.Graph(directed=False)
    >>> create_barabasi_albert(g, n=1000, m=3, seed=42)
    >>> print(g.num_vertices(), g.num_edges())
    1000 2994
    """
    _create_barabasi_albert_raw(g._Graph__graph, n, m, seed)


__all__ = ["create_barabasi_albert"]
