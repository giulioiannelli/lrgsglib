"""Thresholding utilities for weighted graphs."""

from __future__ import annotations

from typing import Tuple

import numpy as np
from networkx import Graph, connected_components

from .base import get_giant_component

__all__ = [
    "compute_threshold_stats",
    "compute_threshold_stats_fast",
    "find_exact_detachment_threshold",
    "select_threshold_and_graph",
    "threshold_graph",
]


def compute_threshold_stats(
    G0: Graph, n_points: int = 0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return threshold values and connectivity statistics for ``G0``."""
    edges_data = [(u, v, data["weight"]) for u, v, data in G0.edges(data=True)]
    weights = np.array([w for _, _, w in edges_data])
    min_weight = weights.min()
    max_weight = weights.max()

    n_points = n_points or G0.number_of_nodes()
    Th = np.logspace(np.log10(min_weight), np.log10(max_weight), n_points)
    Pinf = np.zeros(len(Th))
    Einf = np.zeros(len(Th))

    N0 = G0.number_of_nodes()
    E0 = len(edges_data)
    threshold_indices = np.argsort(Th)[::-1]

    for idx in threshold_indices:
        threshold = Th[idx]
        valid_edges = [(u, v) for u, v, w in edges_data if w >= threshold]

        if not valid_edges:
            Pinf[idx] = 0
            Einf[idx] = 0
            continue

        F = Graph()
        F.add_nodes_from(G0.nodes())
        F.add_edges_from(valid_edges)

        components = list(connected_components(F))
        if not components:
            Pinf[idx] = 0
            Einf[idx] = 0
        else:
            giant_nodes = max(components, key=len)
            giant_size = len(giant_nodes)
            giant_edges = sum(
                1 for u, v in valid_edges if u in giant_nodes and v in giant_nodes
            )
            Pinf[idx] = giant_size / N0
            Einf[idx] = giant_edges / E0

    return Th, Einf, Pinf


def compute_threshold_stats_fast(
    G0: Graph, n_points: int = 0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fast threshold statistics using Union-Find."""
    edges_data = [(u, v, data["weight"]) for u, v, data in G0.edges(data=True)]
    weights = np.array([w for _, _, w in edges_data])

    nodes = list(G0.nodes())
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    n_nodes = len(nodes)

    n_points = n_points or n_nodes
    Th = np.logspace(np.log10(weights.min()), np.log10(weights.max()), n_points)
    Pinf = np.zeros(len(Th))
    Einf = np.zeros(len(Th))

    E0 = len(edges_data)
    edge_indices = np.argsort(weights)[::-1]
    sorted_edges = [
        (edges_data[i][0], edges_data[i][1], weights[i]) for i in edge_indices
    ]

    for i, threshold in enumerate(Th):
        parent = list(range(n_nodes))
        rank = [0] * n_nodes
        component_size = [1] * n_nodes

        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            px, py = find(x), find(y)
            if px == py:
                return
            if rank[px] < rank[py]:
                px, py = py, px
            parent[py] = px
            component_size[px] += component_size[py]
            if rank[px] == rank[py]:
                rank[px] += 1

        valid_edge_count = 0
        for u, v, w in sorted_edges:
            if w >= threshold:
                union(node_to_idx[u], node_to_idx[v])
                valid_edge_count += 1
            else:
                break

        if valid_edge_count == 0:
            Pinf[i] = 0
            Einf[i] = 0
        else:
            max_component_size = max(component_size[find(j)] for j in range(n_nodes))
            giant_root = None
            for j in range(n_nodes):
                if component_size[find(j)] == max_component_size:
                    giant_root = find(j)
                    break

            giant_nodes = {nodes[j] for j in range(n_nodes) if find(j) == giant_root}
            giant_edges = sum(
                1
                for u, v, w in edges_data
                if w >= threshold and u in giant_nodes and v in giant_nodes
            )
            Pinf[i] = max_component_size / n_nodes
            Einf[i] = giant_edges / E0

    return Th, Einf, Pinf


def find_exact_detachment_threshold(corr_mat: np.ndarray) -> float:
    """Return the threshold where the first node detaches from the giant component."""
    n = corr_mat.shape[0]
    triu_indices = np.triu_indices(n, k=1)
    edge_weights = np.abs(corr_mat[triu_indices])
    sorted_weights = np.sort(edge_weights)

    left, right = 0, len(sorted_weights) - 1
    while left < right:
        mid = (left + right) // 2
        threshold = sorted_weights[mid]

        adj_matrix = np.abs(corr_mat) >= threshold
        np.fill_diagonal(adj_matrix, False)

        reach_matrix = adj_matrix.copy().astype(int)
        for _ in range(n - 1):
            reach_matrix = np.dot(reach_matrix, adj_matrix.astype(int))
            if np.any(reach_matrix):
                break

        connected = np.all(reach_matrix + reach_matrix.T + np.eye(n) > 0)

        if connected:
            left = mid + 1
        else:
            right = mid

    return sorted_weights[left] if left < len(sorted_weights) else sorted_weights[-1]


def select_threshold_and_graph(
    G0: Graph, Th: np.ndarray, Pinf: np.ndarray, percentage: float = 0.9
) -> Tuple[float, Graph]:
    """Select optimal threshold and return the thresholded graph."""
    indices = np.where(Pinf < percentage)[0]
    if indices.size > 0:
        best_threshold = Th[indices[0]]
    else:
        best_threshold = Th[-1]

    G_thresh = G0.copy()
    G_thresh.remove_edges_from(
        [(u, v) for u, v, w in G_thresh.edges(data="weight") if w < best_threshold]
    )
    return best_threshold, G_thresh


def threshold_graph(
    G: Graph, percentage: float = 0.9
) -> Tuple[float, Graph, np.ndarray, np.ndarray, np.ndarray]:
    """Return the graph ``G`` thresholded at an optimal value."""
    if any(data["weight"] < 0 for _, _, data in G.edges(data=True)):
        raise ValueError("Graph contains negative weights, which is not allowed.")

    G0 = get_giant_component(G)
    Th, Einf, Pinf = compute_threshold_stats(G0)
    best_threshold, G_thresh = select_threshold_and_graph(G0, Th, Pinf, percentage)

    return best_threshold, G_thresh, Th, Einf, Pinf
