"""Multispectral graph generators and SignedGraph subclass."""

from .MultispectralGraph import MultispectralGraph
from .generators_msg import (
    multiplicative_cascade_probability_matrix,
    multiplicative_cascade_graph,
)

__all__ = [
    "MultispectralGraph",
    "multiplicative_cascade_probability_matrix",
    "multiplicative_cascade_graph",
]
