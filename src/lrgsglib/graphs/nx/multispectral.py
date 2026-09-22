"""Multispectral graph facade module.

Re-exports multispectral graph classes from MultispectralGraphNX,
MultiplicativeCascadeNX, VicsekNX, and HierarchicalModularNX packages
for cleaner imports.
"""

from .HierarchicalModularNX import (
    HierarchicalModularNetworkNX,
)
from .MultiplicativeCascadeNX import (  # Backward compatibility alias
    MultiplicativeCascadeGraph,
    MultiplicativeCascadeGraphNX,
)
from .MultispectralGraphNX import (  # Backward compatibility alias
    MultispectralGraph,
    MultispectralGraphNX,
)
from .VicsekNX import (  # Backward compatibility alias
    VicsekGraph,
    VicsekGraphNX,
)

__all__ = [
    # Base
    "MultispectralGraphNX",
    "MultispectralGraph",
    # Multiplicative Cascade
    "MultiplicativeCascadeGraphNX",
    "MultiplicativeCascadeGraph",
    # Vicsek
    "VicsekGraphNX",
    "VicsekGraph",
    # Hierarchical Modular
    "HierarchicalModularNetworkNX",
]
