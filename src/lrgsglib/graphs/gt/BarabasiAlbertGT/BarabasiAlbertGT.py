"""
BarabasiAlbertGT: graph-tool implementation with C++ backend.

Uses native C++ extension for high-performance graph generation.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
import graph_tool.all as gt

from ..SignedGraphGT import SignedGraphGT
from .cpp import create_barabasi_albert


class BarabasiAlbertGT(SignedGraphGT):
    """
    Barabasi-Albert preferential attachment graph using C++ backend.

    Creates a scale-free network using preferential attachment mechanism.
    New nodes connect to m existing nodes with probability proportional
    to their degree.

    Parameters
    ----------
    n : int
        Number of nodes in final graph.
    m : int
        Number of edges each new node creates (must be <= n).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    n : int
        Number of nodes.
    m : int
        Edges per new node.
    gamma : float
        Expected degree exponent (~3.0 for BA model).

    Examples
    --------
    >>> ba = BarabasiAlbertGT(n=500, m=2, pflip=0.15, seed=42)
    >>> print(f"Nodes: {ba.N}, Edges: {ba.num_edges}")
    >>> hubs = ba.get_hub_nodes(top_k=5)

    Notes
    -----
    The BA model produces scale-free networks with degree distribution
    P(k) ~ k^(-gamma) where gamma ≈ 3.

    Unlike Erdos-Renyi graphs, BA graphs are always connected by construction,
    so no giant component extraction is needed.
    """

    # Expected degree exponent for BA model
    gamma = 3.0

    def __init__(
        self,
        n: int,
        m: int,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ):
        # Validate inputs
        if m > n:
            raise ValueError(f"m ({m}) must be <= n ({n})")
        if m < 1:
            raise ValueError(f"m must be >= 1, got {m}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.m = m

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph using C++ backend
        G = self._generate_graph(seed or 0)

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

        # Apply sign flips if requested
        if pflip > 0:
            self.flip_random_fract_edges()

    def _generate_graph(self, seed: int) -> gt.Graph:
        """Generate Barabasi-Albert graph using C++ extension."""
        G = gt.Graph(directed=False)
        create_barabasi_albert(G, self.n, self.m, seed)

        # Add sign property (all positive initially)
        # Use vectorized initialization for speed
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G

    def get_expected_degree(self) -> float:
        """Return expected average degree."""
        return 2 * self.m

    def get_hub_nodes(self, top_k: int = 10) -> List[Tuple[int, int]]:
        """
        Get the top-k hub nodes by degree.

        Parameters
        ----------
        top_k : int
            Number of top hub nodes to return.

        Returns
        -------
        list of (node_id, degree) tuples
            Sorted by degree descending.
        """
        degrees = [(int(v), v.out_degree()) for v in self.G.vertices()]
        degrees.sort(key=lambda x: x[1], reverse=True)
        return degrees[:top_k]

    def get_degree_distribution(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the degree distribution.

        Returns
        -------
        tuple of (degrees, counts)
            Unique degree values and their frequencies.
        """
        degrees = np.array([v.out_degree() for v in self.G.vertices()])
        unique, counts = np.unique(degrees, return_counts=True)
        return unique, counts

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"N={self.n}_m={self.m}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"BarabasiAlbertGT(n={self.n}, m={self.m}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["BarabasiAlbertGT"]
