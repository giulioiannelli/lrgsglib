"""
Generators for the generalized Sompolinsky-Crisanti-Sommers (SCS) model.
"""

from __future__ import annotations

import math
from typing import Dict, Tuple, Union, Iterable, Optional

import networkx as nx
import numpy as np


EdgeKey = Tuple[int, int]
WeightMap = Dict[EdgeKey, float]
GraphLike = Union[nx.Graph, nx.DiGraph]


def generate_scs_generalized_graph(
    N: int,
    gamma: float,
    J0: float = 0.0,
    J: float = 1.0,
    g: float = 1.0,
    diagonal: Optional[Union[float, Iterable[float]]] = None,
    rng: np.random.Generator | None = None,
) -> GraphLike:
    """
    Build a fully connected graph that follows the generalized SCS ensemble.

    The weight distribution is parameterised by the average coupling strength
    ``J0`` and the coupling variance scale ``J``:

    .. math::
        \\mu = J_0 / N, \\qquad \\sigma^2 = J^2 / N

    Reciprocity between opposing edges is controlled by ``gamma ∈ [-1, 1]``.
    For ``gamma = 1`` the graph is undirected and strictly symmetric, whereas
    ``gamma = -1`` produces perfectly anti-correlated reciprocal weights
    (i.e., ``W_{ji} = 2\\mu - W_{ij}``). Intermediate values yield directed
    graphs where each unordered pair receives Gaussian weights with the desired
    Pearson correlation ``gamma``.

    Parameters
    ----------
    N : int
        Number of nodes in the graph.
    gamma : float
        Correlation parameter in [-1, 1].
    J0 : float, optional
        Mean coupling strength parameter (default: 0.0).
    J : float, optional
        Coupling variance scale parameter (default: 1.0).
    g : float, optional
        Global multiplicative factor applied to every generated weight
        (default: 1.0).
    diagonal : float or Iterable[float], optional
        Controls the diagonal (self-loop) weights. If a scalar is provided,
        all diagonal entries are set to that value. If an iterable of length
        ``N`` is provided, the i-th value is used for node ``i``. If ``None``,
        defaults to legacy behaviour: for ``gamma ≈ 1`` the diagonal is set to
        ``-1``; for other ``gamma`` values, no self-loops are created.
    rng : np.random.Generator, optional
        Random number generator. If None, uses np.random.default_rng().
    """
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    if not -1.0 <= gamma <= 1.0:
        raise ValueError("gamma must lie in the interval [-1, 1].")

    if rng is None:
        rng = np.random.default_rng()

    mu = J0 / N
    sigma = J / math.sqrt(N)

    # Prepare optional diagonal array (None means handled later)
    diag_arr: Optional[np.ndarray]
    if diagonal is None:
        diag_arr = None
    else:
        if isinstance(diagonal, (int, float, np.floating)):
            diag_arr = np.full(N, float(diagonal), dtype=float)
        else:
            diag_arr = np.asarray(list(diagonal), dtype=float)
            if diag_arr.shape != (N,):
                raise ValueError("diagonal must be a scalar or an iterable of length N")

    # Global scaling factor sanity check
    if not isinstance(g, (int, float, np.floating)):
        raise TypeError("g must be a scalar (float or int)")
    g = float(g)

    if np.isclose(gamma, 1.0, atol=1e-12):
        graph: GraphLike = nx.complete_graph(N)
        samples = rng.normal(loc=0.0, scale=1.0, size=(N, N))
        weights: WeightMap = {}
        for i in range(N):
            for j in range(i + 1, N):
                weight = g * (mu + sigma * samples[i, j])
                weights[(i, j)] = weight
        # Assign off-diagonal weights
        nx.set_edge_attributes(graph, weights, "weight")
        # Assign (and create) diagonal self-loops
        for i in range(N):
            if diag_arr is not None:
                dval = float(diag_arr[i])
            else:
                # legacy behaviour for symmetric case
                dval = -1.0
            graph.add_edge(i, i, weight=dval)
        return graph

    graph = nx.complete_graph(N, create_using=nx.DiGraph())
    sqrt_term = math.sqrt(max(0.0, 1.0 - gamma * gamma))
    weights = {}
    for i in range(N):
        for j in range(i + 1, N):
            x = rng.normal()
            y = rng.normal()
            w_ij = g * (mu + sigma * x)
            w_ji = g * (mu + sigma * (gamma * x + sqrt_term * y))
            weights[(i, j)] = w_ij
            weights[(j, i)] = w_ji
    nx.set_edge_attributes(graph, weights, "weight")
    # Optional diagonal (no legacy default for directed case): create self-loops
    if diag_arr is not None:
        for i in range(N):
            graph.add_edge(i, i, weight=float(diag_arr[i]))
    return graph
