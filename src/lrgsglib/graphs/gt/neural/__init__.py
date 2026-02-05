"""
graph-tool neural network graph implementations with GT suffix.

Note: These fall back to NetworkX implementations.
"""

from .HofieldNNGT import HofieldNNGT
from .SCSGeneralizedNNGT import SCSGeneralizedNNGT

__all__ = [
    "HofieldNNGT",
    "SCSGeneralizedNNGT",
]
