"""Neural network graph facade module.

Re-exports neural network graph classes from HofieldNNNX and SCSGeneralizedNNNX
packages for cleaner imports.
"""

from .HofieldNNNX import (
    HofieldNNNX,
    # Backward compatibility alias
    HofieldNN,
    init_mnist_patterns,
)
from .SCSGeneralizedNNNX import SCSGeneralizedNNNX

__all__ = [
    "HofieldNNNX",
    "HofieldNN",
    "init_mnist_patterns",
    "SCSGeneralizedNNNX",
]
