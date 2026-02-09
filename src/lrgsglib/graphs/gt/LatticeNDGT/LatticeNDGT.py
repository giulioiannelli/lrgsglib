"""N-dimensional lattice using graph-tool."""
from __future__ import annotations

try:
    from graph_tool import Graph
    from graph_tool.generation import lattice
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object
    lattice = None


class LatticeNDGT(Graph):
    """N-dimensional lattice graph using graph-tool."""

    def __init__(self, shape: tuple, **kwargs):
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")
        periodic = kwargs.pop('periodic', True)
        g = lattice(shape=shape, periodic=periodic)
        super().__init__(g=g, **kwargs)

__all__ = ["LatticeNDGT"]
