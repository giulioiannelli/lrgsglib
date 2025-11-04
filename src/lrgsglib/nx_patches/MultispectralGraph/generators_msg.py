"""Graph generators for the :mod:`MultispectralGraph` module."""

from __future__ import annotations

import random

import networkx as nx
import numpy as np

from typing import List, Sequence, Tuple


ProbabilityMatrix = np.ndarray
Node2D = Tuple[int, int]


# Rebuild public API after all definitions
__all__ = [
    "multiplicative_cascade_probability_matrix",
    "multiplicative_cascade_graph",
    "multiplicative_cascade_exp_clocks",
    "initial_measure",
    "link_probabilities",
    "palla_lovasz_vicksek_graph",
]


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


def initial_measure(m: int, symmetric: bool = True) -> np.ndarray:
    """Generate a random m x m probability matrix with sum 1.

    If ``symmetric`` is True, the matrix is symmetrized by ``(M + M.T) / 2``
    before normalizing to sum to 1.
    """
    if not isinstance(m, (int, np.integer)) or m <= 0:
        raise ValueError("m must be a positive integer.")
    M = np.random.rand(int(m), int(m)).astype(float)
    if symmetric:
        M = (M + M.T) / 2.0
    s = float(M.sum())
    if s == 0.0:
        raise ValueError("initial_measure produced a zero-sum matrix.")
    M /= s
    return M


def link_probabilities(pij: np.ndarray, k: int) -> np.ndarray:
    """Iterate the generating measure ``pij`` via Kronecker products ``k`` times."""
    if not isinstance(pij, np.ndarray) or pij.ndim != 2 or pij.shape[0] != pij.shape[1]:
        raise ValueError("pij must be a square 2D numpy array.")
    if not isinstance(k, (int, np.integer)) or k <= 0:
        raise ValueError("k must be a positive integer.")
    Pk = pij.copy().astype(float)
    for _ in range(int(k) - 1):
        Pk = np.kron(Pk, pij)
    return Pk


def multiplicative_cascade_probability_matrix(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    iterations: int = 7,
    *,
    stochastic: bool = False,
) -> ProbabilityMatrix:
    """Return the probability matrix obtained via a multiplicative cascade.

    When ``stochastic`` is False, this uses a fixed 2x2 base seed [[p1,p2],[p3,p4]].
    When ``stochastic`` is True, each block placement draws a fresh 2x2 seed
    with entries chosen from {p1,p2,p3,p4} independently.
    """
    if not isinstance(iterations, (int, np.integer)) or iterations < 0:
        raise ValueError("iterations must be a non-negative integer.")

    _validate_probabilities((p1, p2, p3, p4))

    if not stochastic:
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

    # stochastic path
    A = np.random.choice([p1, p2, p3, p4], (2, 2)).astype(float)
    A_new = A.copy()
    for _ in range(iterations):
        A = A_new.copy()
        L_old = A.shape[0]
        L_new = L_old * 2
        A_new = np.zeros((L_new, L_new), dtype=float)
        for i in range(L_old):
            for j in range(L_old):
                seed = np.random.choice([p1, p2, p3, p4], (2, 2)).astype(float)
                A_new[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] += seed * A[i, j]
    return A_new



def multiplicative_cascade_graph(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    *,
    fraction: float,
    iterations: int = 7,
    stochastic: bool = False,
    probabilities: ProbabilityMatrix | None = None,
) -> nx.Graph:
    """Generate a graph sampled from a multiplicative cascade distribution."""
    probabilities = (
        multiplicative_cascade_probability_matrix(
            p1, p2, p3, p4, iterations, stochastic=stochastic
        )
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



def palla_lovasz_vicksek_graph(N: int, pij: np.ndarray, k: int) -> Tuple[nx.Graph, np.ndarray]:
    """Generate a Vicsek graph and return the graph and final probability matrix.

    Steps:
    1) Build the link probability matrix P_k by iterating p_ij k times.
    2) Distribute N nodes uniformly on [0,1].
    3) Connect each pair with probability given by P_k at their coordinates.
    """
    if not isinstance(N, (int, np.integer)) or N <= 0:
        raise ValueError("N must be a positive integer.")

    pij = np.asarray(pij, dtype=float)
    Pk = link_probabilities(pij, k)
    m = pij.shape[0]
    M = int(m ** int(k))

    coords = np.random.rand(int(N))
    indices = np.floor(coords * M).astype(int)
    indices = np.clip(indices, 0, M - 1)

    G = nx.Graph()
    G.add_nodes_from(range(int(N)))

    for i in range(int(N)):
        for j in range(i + 1, int(N)):
            p_link = float(Pk[indices[i], indices[j]])
            if np.random.rand() < p_link:
                G.add_edge(i, j)

    return G, Pk


def multiplicative_cascade_exp_clocks(
    adj: ProbabilityMatrix,
    fraction: float,
) -> nx.Graph:
    """Sample a graph from a precomputed cascade matrix using exponential clocks.

    This implements the selection rule:
    - draw i.i.d. U ~ Uniform(0,1) for each cell
    - compute T = -log(1-U)/A_new (with T = +inf where A_new == 0)
    - select the k cells with smallest T, with k ≈ frazione * N^2
    - induce a subgraph of the periodic 2D grid and return its giant component
    """
    if not isinstance(adj, np.ndarray) or adj.ndim != 2 or adj.shape[0] != adj.shape[1]:
        raise ValueError("A_new must be a square 2D numpy array.")
    if not (0.0 < float(fraction) <= 1.0):
        raise ValueError("frazione must lie in the open interval (0, 1].")

    N = int(adj.shape[0])
    k = int(round(float(fraction) * (N ** 2)))
    k = max(1, min(k, N * N))

    U = np.random.random((N, N))
    # Avoid division by zero: entries with zero probability get infinite time
    with np.errstate(divide="ignore", invalid="ignore"):
        time = -np.log1p(-U) / adj
        time = np.where(adj > 0.0, time, np.inf)

    flat_indices = np.argpartition(time.ravel(), k - 1)[:k]
    indices_2d = np.column_stack(np.unravel_index(flat_indices, (N, N)))
    selected_nodes = [tuple(idx) for idx in indices_2d]

    G = nx.grid_2d_graph(N, N, periodic=True)
    G_sub = G.subgraph(selected_nodes).copy()

    if G_sub.number_of_nodes() == 0:
        raise ValueError("No nodes selected for the cascade graph.")

    if not nx.is_connected(G_sub):
        components = nx.connected_components(G_sub)
        giant_nodes = max(components, key=len)
        return G_sub.subgraph(giant_nodes).copy()
    return G_sub
