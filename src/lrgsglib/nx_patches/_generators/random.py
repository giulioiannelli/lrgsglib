"""
Consolidated random graph generators.

This module provides a unified interface to random graph generation functions,
including both existing types and new implementations.

Existing Types
--------------
- erdos_renyi_graph: Erdos-Renyi random graph
- watts_strogatz_graph: Watts-Strogatz small-world graph

New Types (Phase 3)
-------------------
- barabasi_albert_graph: Scale-free network
- stochastic_block_model: Community-structured network
- k_regular_graph: Random k-regular graph
- configuration_model: Graph with prescribed degree sequence
- random_geometric_graph: Spatial random graph

Example
-------
>>> from lrgsglib.nx_patches._generators import random
>>> G = random.erdos_renyi_graph(100, 0.1)
>>> G.number_of_nodes()
100
"""

import networkx as nx
from networkx import Graph
from typing import Optional, Sequence
import numpy as np

__all__ = [
    # Core random graph generators
    "erdos_renyi_graph",
    "watts_strogatz_graph",
    # Scale-free and community models
    "barabasi_albert_graph",
    "stochastic_block_model",
    # Other random models
    "k_regular_graph",
    "configuration_model",
    "random_geometric_graph",
    # Utilities
    "extract_giant_component",
]


# =============================================================================
# Core Random Graph Generators
# =============================================================================

def erdos_renyi_graph(
    n: int,
    p: float,
    seed: Optional[int] = None,
    directed: bool = False,
) -> Graph:
    """
    Generate an Erdos-Renyi random graph G(n, p).

    Parameters
    ----------
    n : int
        Number of nodes.
    p : float
        Probability of edge creation.
    seed : int, optional
        Random seed for reproducibility.
    directed : bool, default False
        Whether to create a directed graph.

    Returns
    -------
    Graph
        A NetworkX graph.

    Example
    -------
    >>> G = erdos_renyi_graph(100, 0.05, seed=42)
    >>> G.number_of_nodes()
    100
    """
    return nx.erdos_renyi_graph(n, p, seed=seed, directed=directed)


def watts_strogatz_graph(
    n: int,
    k: int,
    p: float,
    seed: Optional[int] = None,
) -> Graph:
    """
    Generate a Watts-Strogatz small-world graph.

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Each node is connected to k nearest neighbors in ring topology.
    p : float
        Probability of rewiring each edge.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Graph
        A NetworkX graph with small-world properties.

    Example
    -------
    >>> G = watts_strogatz_graph(100, 4, 0.3, seed=42)
    >>> G.number_of_nodes()
    100
    """
    return nx.watts_strogatz_graph(n, k, p, seed=seed)


# =============================================================================
# Scale-Free and Community Models
# =============================================================================

def barabasi_albert_graph(
    n: int,
    m: int,
    seed: Optional[int] = None,
) -> Graph:
    """
    Generate a Barabasi-Albert scale-free network.

    Parameters
    ----------
    n : int
        Number of nodes.
    m : int
        Number of edges to attach from a new node to existing nodes.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Graph
        A NetworkX graph with power-law degree distribution.

    Notes
    -----
    Scale-free networks exhibit different localization properties
    compared to random graphs due to the presence of hubs.

    Example
    -------
    >>> G = barabasi_albert_graph(100, 2, seed=42)
    >>> G.number_of_nodes()
    100
    """
    return nx.barabasi_albert_graph(n, m, seed=seed)


def stochastic_block_model(
    sizes: Sequence[int],
    p_matrix: Sequence[Sequence[float]],
    seed: Optional[int] = None,
) -> Graph:
    """
    Generate a stochastic block model graph.

    Parameters
    ----------
    sizes : sequence of int
        Number of nodes in each block/community.
    p_matrix : sequence of sequence of float
        Probability matrix where p_matrix[i][j] is the probability
        of an edge between nodes in block i and block j.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Graph
        A NetworkX graph with community structure.

    Notes
    -----
    SBM is useful for studying spectral gaps in modular networks
    and community detection algorithms.

    Example
    -------
    >>> sizes = [25, 25, 25, 25]
    >>> p_in, p_out = 0.3, 0.05
    >>> p_matrix = [[p_in if i == j else p_out for j in range(4)] for i in range(4)]
    >>> G = stochastic_block_model(sizes, p_matrix, seed=42)
    >>> G.number_of_nodes()
    100
    """
    return nx.stochastic_block_model(sizes, p_matrix, seed=seed)


# =============================================================================
# Other Random Models
# =============================================================================

def k_regular_graph(
    n: int,
    k: int,
    seed: Optional[int] = None,
) -> Graph:
    """
    Generate a random k-regular graph.

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Degree of each node.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Graph
        A NetworkX graph where every node has degree k.

    Notes
    -----
    Regular graphs serve as baseline for studying disorder effects
    since degree heterogeneity is eliminated.

    Example
    -------
    >>> G = k_regular_graph(100, 4, seed=42)
    >>> all(d == 4 for _, d in G.degree())
    True
    """
    return nx.random_regular_graph(k, n, seed=seed)


def configuration_model(
    degree_sequence: Sequence[int],
    seed: Optional[int] = None,
    create_using: Optional[Graph] = None,
) -> Graph:
    """
    Generate a random graph with prescribed degree sequence.

    Parameters
    ----------
    degree_sequence : sequence of int
        The degree sequence to match.
    seed : int, optional
        Random seed for reproducibility.
    create_using : Graph, optional
        Graph type to create.

    Returns
    -------
    Graph
        A NetworkX graph with the given degree sequence.

    Notes
    -----
    Useful for decoupling degree distribution effects from other
    structural properties in spectral analysis.

    Example
    -------
    >>> degrees = [3, 3, 3, 3, 2, 2]
    >>> G = configuration_model(degrees, seed=42)
    """
    G = nx.configuration_model(degree_sequence, seed=seed, create_using=create_using)
    # Remove self-loops and multi-edges
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    return G


def random_geometric_graph(
    n: int,
    radius: float,
    dim: int = 2,
    seed: Optional[int] = None,
) -> Graph:
    """
    Generate a random geometric graph.

    Parameters
    ----------
    n : int
        Number of nodes.
    radius : float
        Distance threshold for edge creation.
    dim : int, default 2
        Dimension of the embedding space.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Graph
        A NetworkX graph where nodes are connected if within radius.

    Notes
    -----
    Random geometric graphs exhibit spatial constraints and are
    relevant for percolation studies.

    Example
    -------
    >>> G = random_geometric_graph(100, 0.2, dim=2, seed=42)
    """
    return nx.random_geometric_graph(n, radius, dim=dim, seed=seed)


# =============================================================================
# Utilities
# =============================================================================

def extract_giant_component(G: Graph) -> Graph:
    """
    Extract the largest connected component from a graph.

    Parameters
    ----------
    G : Graph
        Input graph.

    Returns
    -------
    Graph
        The largest connected component as a new graph.

    Example
    -------
    >>> G = erdos_renyi_graph(100, 0.01, seed=42)
    >>> GC = extract_giant_component(G)
    >>> nx.is_connected(GC)
    True
    """
    if nx.is_connected(G):
        return G.copy()
    components = nx.connected_components(G)
    largest = max(components, key=len)
    return G.subgraph(largest).copy()
