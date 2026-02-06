"""
Temporal graph module.

This module provides classes for time-varying/dynamic networks:
- TemporalGraph: Base class for temporal networks with snapshots
- TemporalSignedGraph: Temporal networks with signed edges

Example
-------
>>> from lrgsglib.nx_patches.temporal import TemporalGraph
>>> tg = TemporalGraph(n_nodes=100, n_snapshots=10)
>>> tg.add_snapshot(0, edges=[(0, 1), (1, 2)])
>>> G_t5 = tg.get_snapshot(5)
"""

from .temporal_graph import TemporalGraph, TemporalSignedGraph

__all__ = [
    "TemporalGraph",
    "TemporalSignedGraph",
]
