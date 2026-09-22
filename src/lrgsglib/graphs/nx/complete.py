"""Complete graph facade module.

Re-exports complete graph classes from CompleteGraphNX and FullyConnectedNX
packages for cleaner imports.
"""

from .CompleteGraphNX import (  # Backward compatibility alias
    CompleteGraph,
    CompleteGraphNX,
)
from .FullyConnectedNX import (  # Backward compatibility alias
    FullyConnected,
    FullyConnectedNX,
)

__all__ = [
    "CompleteGraphNX",
    "CompleteGraph",
    "FullyConnectedNX",
    "FullyConnected",
]
