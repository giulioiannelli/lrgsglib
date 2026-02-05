"""Helpers to query neighbours and cycles in graphs."""

from __future__ import annotations

import warnings
from collections import deque
from typing import Any, List

import numpy as np
from networkx import (
    Graph,
    NetworkXError,
    is_connected,
    single_source_shortest_path_length,
)

__all__ = [
    "get_kth_order_neighbours",
    "get_neighbors_at_distance",
    "get_smallest_cycle_graph_node",
]


def get_kth_order_neighbours(
    G: Graph, node: Any, order: int = 1
) -> List[Any]:
    """Return nodes exactly ``order`` steps away from ``node``."""
    if not isinstance(order, int) or order <= 0:
        raise ValueError("The `order` must be a positive integer.")

    if node not in G:
        raise NetworkXError(f"The node {node} is not in the graph.")

    if not is_connected(G):
        warnings.warn(
            "The graph is disconnected. Results might be meaningless.",
            UserWarning,
        )

    path_lengths = single_source_shortest_path_length(G, node, cutoff=order)
    return [n for n, d in path_lengths.items() if d == order]


def get_neighbors_at_distance(G: Graph, node: Any, distance: int) -> List[Any]:
    """Return neighbors of ``node`` at given ``distance``."""
    n_dict = single_source_shortest_path_length(G, node)
    return [n for n, d in n_dict.items() if d == distance]


def get_smallest_cycle_graph_node(G: Graph, start_node: Any) -> List[Any]:
    """Return the smallest cycle containing ``start_node`` if any."""
    visited = {start_node: None}
    queue = deque([(start_node, None)])
    smallest_cycle = None

    while queue:
        current_node, parent = queue.popleft()

        for neighbor in G.neighbors(current_node):
            if neighbor == parent:
                continue
            if neighbor in visited:
                cycle = []
                node = current_node
                while node is not None and node != neighbor:
                    cycle.append(node)
                    node = visited[node]
                cycle.append(neighbor)
                cycle.append(current_node)
                if smallest_cycle is None or len(cycle) < len(smallest_cycle):
                    smallest_cycle = cycle
                    if len(smallest_cycle) == 3:
                        return smallest_cycle
            else:
                visited[neighbor] = current_node
                queue.append((neighbor, current_node))

    return smallest_cycle
