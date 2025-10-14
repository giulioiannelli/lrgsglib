"""Graph generators for the :mod:`MultispectralGraph` module."""

from __future__ import annotations

from typing import List, Sequence, Tuple

import networkx as nx
import numpy as np

from ..common import random

ProbabilityMatrix = np.ndarray
Node2D = Tuple[int, int]


def _validate_probabilities(probabilities: Sequence[float]) -> None:
    """Validate that probabilities are numeric and lie in ``[0, 1]``."""
    if len(probabilities) != 4:
        raise ValueError("Exactly four probabilities are required for the cascade seed.")
    for idx, value in enumerate(probabilities):
        if not isinstance(value, (int, float, np.floating)):
            raise TypeError(
                f"Probability p{idx + 1} must be a real number, received {type(value)!r}."
            )
        if not (0.0 <= float(value) <= 1.0):
            raise ValueError(
                f"Probability p{idx + 1} must be within [0, 1], received {value}."
            )


def multiplicative_cascade_probability_matrix(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    iterations: int = 7,
) -> ProbabilityMatrix:
    """Return the probability matrix obtained via a multiplicative cascade."""
    if not isinstance(iterations, (int, np.integer)) or iterations < 0:
        raise ValueError("iterations must be a non-negative integer.")

    _validate_probabilities((p1, p2, p3, p4))

    base = np.array([[p1, p2], [p3, p4]], dtype=float)
    matrix = base.copy()

    for _ in range(iterations):
        previous = matrix
        size = previous.shape[0]
        expanded = np.zeros((size * 2, size * 2), dtype=float)
        for i in range(size):
            for j in range(size):
                expanded[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] += base * previous[i, j]
        matrix = expanded

    return matrix


def _select_nodes(
    graph: nx.Graph,
    probabilities: ProbabilityMatrix,
    fraction: float,
) -> List[Node2D]:
    """Randomly select nodes according to ``probabilities`` until ``fraction`` is met."""
    if not (0.0 < fraction <= 1.0):
        raise ValueError("fraction must lie in the open interval (0, 1].")

    nodes: List[Node2D] = list(graph.nodes())
    total_nodes = len(nodes)
    target = max(1, int(fraction * total_nodes))
    selected: List[Node2D] = []
    attempts = 0
    max_attempts = 10 * total_nodes if total_nodes else 0

    while len(selected) < target and nodes:
        node = random.choice(nodes)
        i, j = node
        if np.random.random() < probabilities[i, j]:
            selected.append(node)
            nodes.remove(node)
        attempts += 1
        if attempts > max_attempts:
            break

    if len(selected) < target and nodes:
        remaining = sorted(nodes, key=lambda node: probabilities[node], reverse=True)
        for node in remaining:
            selected.append(node)
            if len(selected) >= target:
                break

    if not selected:
        raise ValueError("Unable to select any node from the cascade probabilities.")

    return selected


def multiplicative_cascade_graph(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    *,
    fraction: float,
    iterations: int = 7,
    probabilities: ProbabilityMatrix | None = None,
) -> nx.Graph:
    """Generate a graph sampled from a multiplicative cascade distribution."""
    probabilities = (
        multiplicative_cascade_probability_matrix(p1, p2, p3, p4, iterations)
        if probabilities is None
        else probabilities
    )
    size = probabilities.shape[0]

    base_graph = nx.grid_2d_graph(size, size)
    selected_nodes = _select_nodes(base_graph, probabilities, fraction)
    subgraph = base_graph.subgraph(selected_nodes).copy()

    if subgraph.number_of_nodes() == 0:
        raise ValueError("The multiplicative cascade graph is empty.")

    if not nx.is_connected(subgraph):
        giant_component_nodes = max(nx.connected_components(subgraph), key=len)
        subgraph = subgraph.subgraph(giant_component_nodes).copy()

    return subgraph

__all__ = [
    "multiplicative_cascade_probability_matrix",
    "multiplicative_cascade_graph",
]
