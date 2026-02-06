"""Bipartite graph facade module.

Re-exports bipartite graph classes from BipartiteGraphNX and
BipartiteFromDegreeSequenceNX packages for cleaner imports.
"""

from .BipartiteGraphNX import (
    BipartiteGraphNX,
    # Backward compatibility alias
    BipartiteGraph,
)
from .BipartiteFromDegreeSequenceNX import (
    BipartiteFromDegreeSequenceNX,
    # Backward compatibility alias
    BipartiteFromDegreeSequence,
)

__all__ = [
    "BipartiteGraphNX",
    "BipartiteGraph",
    "BipartiteFromDegreeSequenceNX",
    "BipartiteFromDegreeSequence",
]
