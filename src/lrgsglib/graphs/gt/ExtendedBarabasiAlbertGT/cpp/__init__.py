"""
C++ extension for ExtendedBarabasiAlbertGT graph generation.

Provides high-performance Extended Barabasi-Albert preferential attachment
graph construction with initial attractiveness parameter.

Build with: make cpp-make (from lrgsglib root)
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import graph_tool

# Import the compiled extension
_SO_PATH = Path(__file__).parent / "extended_barabasi_albert.so"

try:
    from . import extended_barabasi_albert as _eba_module
    _create_extended_barabasi_albert_raw = _eba_module.create_extended_barabasi_albert
except ImportError as e:
    _eba_module = None
    _import_error = e

    def _create_extended_barabasi_albert_raw(*args, **kwargs):
        raise ImportError(
            f"extended_barabasi_albert C++ extension not available at {_SO_PATH}. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from _import_error


def create_extended_barabasi_albert(
    g: "graph_tool.Graph",
    n: int,
    m: int,
    a: float = 0.0,
    seed: int = 0,
) -> None:
    """
    Create an Extended Barabasi-Albert graph with initial attractiveness.

    Attachment probability: P(i) ~ (k_i + a)
    When a=0, reduces to standard BA. Higher a leads to more uniform degrees.

    Parameters
    ----------
    g : graph_tool.Graph
        An empty graph-tool Graph object to populate.
    n : int
        Total number of nodes.
    m : int
        Number of edges each new node creates.
    a : float, optional
        Initial attractiveness parameter. Default 0.0 (standard BA).
    seed : int, optional
        Random seed for reproducibility. 0 = use random device.

    Examples
    --------
    >>> import graph_tool.all as gt
    >>> from lrgsglib.graphs.gt.ExtendedBarabasiAlbertGT.cpp import create_extended_barabasi_albert
    >>> g = gt.Graph(directed=False)
    >>> create_extended_barabasi_albert(g, n=1000, m=3, a=1.0, seed=42)
    >>> print(g.num_vertices(), g.num_edges())
    1000 2994
    """
    _create_extended_barabasi_albert_raw(g._Graph__graph, n, m, a, seed)


__all__ = ["create_extended_barabasi_albert"]
