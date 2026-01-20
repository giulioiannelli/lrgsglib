"""Hierarchical Modular Network generators.

Implements the hierarchical modular network model where nodes are organized
into nested modules with decreasing connection probabilities at higher
hierarchical levels.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
from typing import Optional


def hierarchical_modular_network(
    levels: int = 3,
    branching: int = 2,
    leaf_nodes: int = 8,
    p_intra: float = 0.8,
    p_ratio: float = 0.1,
    seed: Optional[int] = None,
) -> nx.Graph:
    """Generate a hierarchical modular network.

    Creates a graph with nested community structure. At the lowest level,
    nodes form tight communities with connection probability p_intra.
    At each higher level, communities merge into larger modules with
    connection probability scaled by p_ratio^level.

    Parameters
    ----------
    levels : int
        Number of hierarchical levels (depth of tree). Default: 3
    branching : int
        Branching factor at each level (number of submodules per module).
        Default: 2
    leaf_nodes : int
        Number of nodes in each leaf module (lowest level). Default: 8
    p_intra : float
        Edge probability within leaf modules. Default: 0.8
    p_ratio : float
        Ratio for edge probability decay between levels.
        p(level k) = p_intra * p_ratio^k. Default: 0.1
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    nx.Graph
        Hierarchical modular network with integer node labels.

    Examples
    --------
    >>> G = hierarchical_modular_network(levels=2, branching=2, leaf_nodes=4)
    >>> G.number_of_nodes()
    16
    >>> # 2^2 = 4 leaf modules, each with 4 nodes = 16 total

    Notes
    -----
    Total nodes = branching^levels * leaf_nodes

    The hierarchy is structured as a tree:
    - Root level (k=levels): entire graph
    - Intermediate levels: nested modules
    - Leaf level (k=0): tightest communities

    Connection probability between nodes depends on their hierarchical
    distance: p(edge) = p_intra * p_ratio^k where k is the number of
    levels up to their lowest common ancestor.
    """
    rng = np.random.default_rng(seed)

    n_modules = branching ** levels
    n_nodes = n_modules * leaf_nodes

    G = nx.Graph()
    G.add_nodes_from(range(n_nodes))

    def get_leaf_module(node_id: int) -> int:
        """Get leaf module index for a node."""
        return node_id // leaf_nodes

    def hierarchical_distance(node_i: int, node_j: int) -> int:
        """Return hierarchical distance (0=same leaf module, levels=different root)."""
        mod_i = get_leaf_module(node_i)
        mod_j = get_leaf_module(node_j)

        if mod_i == mod_j:
            return 0  # Same leaf module

        # Find first level where they share a parent
        for lvl in range(1, levels + 1):
            # At level lvl, modules are grouped in blocks of branching^lvl
            block_size = branching ** lvl
            if mod_i // block_size == mod_j // block_size:
                return lvl

        return levels  # Different at all levels

    # Add edges based on hierarchical distance
    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            dist = hierarchical_distance(i, j)
            p_edge = p_intra * (p_ratio ** dist)
            if rng.random() < p_edge:
                G.add_edge(i, j, weight=1.0)

    return G
