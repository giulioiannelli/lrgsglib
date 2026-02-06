"""GraphOfGraphsGT: graph-tool-based hierarchical graph structure.

This module provides the GraphOfGraphsGT class for creating generalized
"graph of graphs" structures where arbitrary graph types can be used
as base and fiber components.
"""

from .GraphOfGraphsGT import GraphOfGraphsGT, GraphOfGraphs
from ._policies import AnchorPolicy, get_anchor_index, resolve_anchor_indices

__all__ = [
    "GraphOfGraphsGT",
    "GraphOfGraphs",
    "AnchorPolicy",
    "get_anchor_index",
    "resolve_anchor_indices",
]
