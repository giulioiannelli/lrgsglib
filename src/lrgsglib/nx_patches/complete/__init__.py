"""
Complete graph family module.

This module provides a clean class hierarchy for fully connected graphs:
- CompleteGraph: Abstract base class for complete graphs
- FullyConnected: Standard fully connected graph

Example
-------
>>> from lrgsglib.nx_patches.complete import FullyConnected
>>> fc = FullyConnected(N=50, pflip=0.1)
>>> fc.flip_random_fract_edges()
"""

from .complete_graph import CompleteGraph
from .fully_connected import FullyConnected

__all__ = [
    "CompleteGraph",
    "FullyConnected",
]
