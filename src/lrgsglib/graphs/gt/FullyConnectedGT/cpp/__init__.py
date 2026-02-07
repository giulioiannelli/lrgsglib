"""
C++ extension for FullyConnectedGT graph generation.

Provides high-performance complete graph construction using graph-tool's C++ backend.

Build with: make cpp-make (from lrgsglib root)
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import graph_tool

# Import the compiled extension
_SO_PATH = Path(__file__).parent / "complete_graph.so"

try:
    from . import complete_graph as _cg_module
    _create_complete_graph_raw = _cg_module.create_complete_graph
except ImportError as e:
    _cg_module = None
    _import_error = e

    def _create_complete_graph_raw(*args, **kwargs):
        raise ImportError(
            f"complete_graph C++ extension not available at {_SO_PATH}. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from _import_error


def create_complete_graph(
    g: "graph_tool.Graph",
    n: int,
) -> None:
    """
    Create a complete (fully connected) graph.

    Every node is connected to every other node.
    Total edges: n*(n-1)/2

    Parameters
    ----------
    g : graph_tool.Graph
        An empty graph-tool Graph object to populate.
    n : int
        Number of nodes.

    Examples
    --------
    >>> import graph_tool.all as gt
    >>> from lrgsglib.graphs.gt.FullyConnectedGT.cpp import create_complete_graph
    >>> g = gt.Graph(directed=False)
    >>> create_complete_graph(g, n=100)
    >>> print(g.num_vertices(), g.num_edges())
    100 4950
    """
    _create_complete_graph_raw(g._Graph__graph, n)


__all__ = ["create_complete_graph"]
