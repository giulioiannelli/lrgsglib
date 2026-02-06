"""Weighted graph wrapper for graph-tool."""

try:
    from graph_tool import Graph
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object


class WeightedGraph(Graph):
    """Weighted graph with edge weight property."""

    def __init__(self, *args, **kwargs):
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")
        super().__init__(*args, **kwargs)
