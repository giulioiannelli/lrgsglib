"""
C++ extension for HolmeKimGT graph generation.

Provides high-performance Holme-Kim powerlaw cluster graph
generation using graph-tool's C++ backend.

Build with: make cpp-make (from lrgsglib root)
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import graph_tool

# Import the compiled extension
_SO_PATH = Path(__file__).parent / "holme_kim.so"

try:
    from . import holme_kim as _hk_module
    _create_holme_kim_raw = _hk_module.create_holme_kim
except ImportError as e:
    _hk_module = None
    _import_error = e

    def _create_holme_kim_raw(*args, **kwargs):
        raise ImportError(
            f"holme_kim C++ extension not available at {_SO_PATH}. "
            "Build with 'make cpp-make' from the lrgsglib root directory."
        ) from _import_error


def create_holme_kim(
    g: "graph_tool.Graph",
    n: int,
    m: int,
    p: float,
    seed: int = 0,
) -> None:
    """
    Create a Holme-Kim powerlaw cluster graph.

    This produces scale-free networks with tunable clustering coefficient.
    It's a BA variant with an additional triangle formation step.

    Parameters
    ----------
    g : graph_tool.Graph
        An empty graph-tool Graph object to populate.
    n : int
        Total number of nodes.
    m : int
        Number of edges each new node creates (must be < n).
    p : float
        Probability of triad formation step (0=pure BA, 1=max clustering).
    seed : int, optional
        Random seed for reproducibility. 0 = use random device.

    Examples
    --------
    >>> import graph_tool.all as gt
    >>> from lrgsglib.graphs.gt.HolmeKimGT.cpp import create_holme_kim
    >>> g = gt.Graph(directed=False)
    >>> create_holme_kim(g, n=1000, m=3, p=0.8, seed=42)
    >>> # This graph will have higher clustering than pure BA
    """
    _create_holme_kim_raw(g._Graph__graph, n, m, p, seed)


__all__ = ["create_holme_kim"]
