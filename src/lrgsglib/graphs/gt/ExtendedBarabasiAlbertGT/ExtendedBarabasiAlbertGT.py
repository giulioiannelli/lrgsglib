"""
ExtendedBarabasiAlbertGT: graph-tool implementation with C++ backend.

Uses native C++ extension for high-performance graph generation.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
import graph_tool.all as gt

from ....gt_patches.signed_graph_gt import SignedGraphGT
from .cpp import create_extended_barabasi_albert


class ExtendedBarabasiAlbertGT(SignedGraphGT):
    """
    Extended Barabasi-Albert preferential attachment graph using C++ backend.

    Attachment probability: P(i) ~ (k_i + a)
    When a=0, this reduces to standard BA.
    Higher a leads to more uniform degree distribution.

    Parameters
    ----------
    n : int
        Number of nodes in final graph.
    m : int
        Number of edges each new node creates.
    a : float
        Initial attractiveness parameter (>= 0).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Examples
    --------
    >>> eba = ExtendedBarabasiAlbertGT(n=500, m=2, a=1.0, pflip=0.15, seed=42)
    >>> eba.flip_random_fract_edges()
    >>> print(f"Nodes: {eba.N}, Edges: {eba.num_edges}")
    """

    def __init__(
        self,
        n: int,
        m: int,
        a: float = 0.0,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ):
        # Validate inputs
        if m > n:
            raise ValueError(f"m ({m}) must be <= n ({n})")
        if m < 1:
            raise ValueError(f"m must be >= 1, got {m}")
        if a < 0:
            raise ValueError(f"a must be >= 0, got {a}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.m = m
        self.a = a

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph using C++ backend
        G = self._generate_graph(seed or 0)

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

    def _generate_graph(self, seed: int) -> gt.Graph:
        """Generate Extended BA graph using C++ extension."""
        G = gt.Graph(directed=False)
        create_extended_barabasi_albert(G, self.n, self.m, self.a, seed)

        # Add sign property (all positive initially)
        # Use vectorized initialization for speed
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G

    def get_expected_degree(self) -> float:
        """Return expected average degree."""
        return 2 * self.m

    def get_hub_nodes(self, top_k: int = 10) -> List[Tuple[int, int]]:
        """Get the top-k hub nodes by degree."""
        degrees = [(int(v), v.out_degree()) for v in self.G.vertices()]
        degrees.sort(key=lambda x: x[1], reverse=True)
        return degrees[:top_k]

    def get_degree_distribution(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get the degree distribution."""
        degrees = np.array([v.out_degree() for v in self.G.vertices()])
        unique, counts = np.unique(degrees, return_counts=True)
        return unique, counts

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"N={self.n}_m={self.m}_a={self.a}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"ExtendedBarabasiAlbertGT(n={self.n}, m={self.m}, a={self.a}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["ExtendedBarabasiAlbertGT"]
