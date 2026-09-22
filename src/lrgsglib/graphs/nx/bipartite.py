"""Bipartite graph facade module.

Re-exports bipartite graph classes from BipartiteGraphNX and
BipartiteFromDegreeSequenceNX packages for cleaner imports.
"""

from .BipartiteFromDegreeSequenceNX import (  # Backward compatibility alias
    BipartiteFromDegreeSequence,
    BipartiteFromDegreeSequenceNX,
)
from .BipartiteGraphNX import (  # Backward compatibility alias
    BipartiteGraph,
    BipartiteGraphNX,
)

__all__ = [
    "BipartiteGraphNX",
    "BipartiteGraph",
    "BipartiteFromDegreeSequenceNX",
    "BipartiteFromDegreeSequence",
]
