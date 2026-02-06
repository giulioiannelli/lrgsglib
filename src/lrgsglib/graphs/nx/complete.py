"""Complete graph facade module.

Re-exports complete graph classes from CompleteGraphNX and FullyConnectedNX
packages for cleaner imports.
"""

from .CompleteGraphNX import (
    CompleteGraphNX,
    # Backward compatibility alias
    CompleteGraph,
)
from .FullyConnectedNX import (
    FullyConnectedNX,
    # Backward compatibility alias
    FullyConnected,
)

__all__ = [
    "CompleteGraphNX",
    "CompleteGraph",
    "FullyConnectedNX",
    "FullyConnected",
]
