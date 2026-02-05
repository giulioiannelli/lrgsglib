"""Multispectral graph facade module.

Re-exports multispectral graph classes from MultispectralGraphNX,
MultiplicativeCascadeNX, VicsekNX, and HierarchicalModularNX packages
for cleaner imports.
"""

from .MultispectralGraphNX import (
    MultispectralGraphNX,
    # Backward compatibility alias
    MultispectralGraph,
)
from .MultiplicativeCascadeNX import (
    MultiplicativeCascadeGraphNX,
    # Backward compatibility alias
    MultiplicativeCascadeGraph,
)
from .VicsekNX import (
    VicsekGraphNX,
    # Backward compatibility alias
    VicsekGraph,
)
from .HierarchicalModularNX import (
    HierarchicalModularNetworkNX,
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
