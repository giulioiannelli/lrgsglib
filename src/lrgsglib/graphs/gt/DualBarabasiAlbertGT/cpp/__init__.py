"""
C++ extension for DualBarabasiAlbertGT graph generation.

Provides high-performance Dual Barabasi-Albert graph construction
with two attachment modes.

Build with: make cpp-make (from lrgsglib root)
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import graph_tool

# Import the compiled extension
_SO_PATH = Path(__file__).parent / "dual_barabasi_albert.so"

try:
    from . import dual_barabasi_albert as _dba_module
    _create_dual_barabasi_albert_raw = _dba_module.create_dual_barabasi_albert
except ImportError as e:
    _dba_module = None
    _import_error = e

    def _create_dual_barabasi_albert_raw(*args, **kwargs):
        raise ImportError(
            f"dual_barabasi_albert C++ extension not available at {_SO_PATH}. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from _import_error


def create_dual_barabasi_albert(
    g: "graph_tool.Graph",
    n: int,
    m1: int,
    m2: int,
    p: float = 0.5,
    seed: int = 0,
) -> None:
    """
    Create a Dual Barabasi-Albert graph with two attachment modes.

    With probability p, new nodes use m1 edges; otherwise m2 edges.

    Parameters
    ----------
    g : graph_tool.Graph
        An empty graph-tool Graph object to populate.
    n : int
        Total number of nodes.
    m1 : int
        Number of edges in mode 1.
    m2 : int
        Number of edges in mode 2.
    p : float, optional
        Probability of using m1 vs m2. Default 0.5.
    seed : int, optional
        Random seed for reproducibility. 0 = use random device.

    Examples
    --------
    >>> import graph_tool.all as gt
    >>> from lrgsglib.graphs.gt.DualBarabasiAlbertGT.cpp import create_dual_barabasi_albert
    >>> g = gt.Graph(directed=False)
    >>> create_dual_barabasi_albert(g, n=1000, m1=2, m2=5, p=0.7, seed=42)
    >>> print(g.num_vertices())
    1000
    """
    _create_dual_barabasi_albert_raw(g._Graph__graph, n, m1, m2, p, seed)


__all__ = ["create_dual_barabasi_albert"]
