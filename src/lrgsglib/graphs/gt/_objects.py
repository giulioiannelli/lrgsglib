"""
gt objects module.

Common imports for graph-tool graphs.
"""
try:
    from graph_tool.all import Graph, lattice
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = None
    lattice = None
