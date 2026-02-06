"""GraphOfGraphsNX: NetworkX-based hierarchical graph structure.

This module provides the GraphOfGraphsNX class for creating generalized
"graph of graphs" structures where arbitrary graph types can be used
as base and fiber components.
"""

from .GraphOfGraphsNX import GraphOfGraphsNX, GraphOfGraphs
from ._policies import AnchorPolicy, get_anchor_index, resolve_anchor_indices

__all__ = [
    "GraphOfGraphsNX",
    "GraphOfGraphs",
    "AnchorPolicy",
    "get_anchor_index",
    "resolve_anchor_indices",
]
