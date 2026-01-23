"""
Fractal graph family module.

This module provides a clean class hierarchy for fractal/self-similar graphs:
- FractalGraph: Abstract base class for fractal graphs
- DGMgraph: Dorogovtsev-Goltsev-Mendes graph
- SierpinskiGraph: Sierpinski gasket/triangle/carpet/tetrahedron

Example
-------
>>> from lrgsglib.nx_patches.fractal import DGMgraph, SierpinskiGraph
>>> dgm = DGMgraph(n=5, pflip=0.1)
>>> dgm.flip_random_fract_edges()
>>> sierp = SierpinskiGraph(n=4, variant="gasket")
"""

from .fractal_graph import FractalGraph
from .dgm_graph import DGMgraph
from .sierpinski import SierpinskiGraph

__all__ = [
    "FractalGraph",
    "DGMgraph",
    "SierpinskiGraph",
]
