"""
graph-tool fractal graph implementations with GT suffix.

Note: These fall back to NetworkX implementations.
"""

from .FractalGraphGT import FractalGraphGT
from .DGMgraphGT import DGMgraphGT
from .SierpinskiGraphGT import SierpinskiGraphGT

__all__ = [
    "FractalGraphGT",
    "DGMgraphGT",
    "SierpinskiGraphGT",
]
