"""Neural network graph facade module.

Re-exports neural network graph classes from HofieldNNNX and SCSGeneralizedNNNX
packages for cleaner imports.
"""

from .HofieldNNNX import (  # Backward compatibility alias
    HofieldNN,
    HofieldNNNX,
    init_mnist_patterns,
)
from .SCSGeneralizedNNNX import SCSGeneralizedNNNX

__all__ = [
    "HofieldNNNX",
    "HofieldNN",
    "init_mnist_patterns",
    "SCSGeneralizedNNNX",
]
