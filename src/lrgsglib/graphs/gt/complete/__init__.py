"""
graph-tool complete graph implementations with GT suffix.

Note: These fall back to NetworkX implementations.
"""

from .CompleteGraphGT import CompleteGraphGT
from .FullyConnectedGT import FullyConnectedGT

__all__ = [
    "CompleteGraphGT",
    "FullyConnectedGT",
]
