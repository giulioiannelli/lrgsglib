"""
DualBarabasiAlbertGT: graph-tool implementation with C++ backend.

Uses native C++ extension for high-performance graph generation.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
import graph_tool.all as gt

from ..SignedGraphGT import SignedGraphGT
from .cpp import create_dual_barabasi_albert


class DualBarabasiAlbertGT(SignedGraphGT):
    """
    Dual Barabasi-Albert graph using C++ backend.

    With probability p, new nodes use m1 edges.
    With probability 1-p, new nodes use m2 edges.

    Parameters
    ----------
    n : int
        Number of nodes in final graph.
    m1 : int
        Number of edges in mode 1.
    m2 : int
        Number of edges in mode 2.
    p : float
        Probability of using m1 (vs m2).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Examples
    --------
    >>> dba = DualBarabasiAlbertGT(n=500, m1=2, m2=5, p=0.7, pflip=0.15, seed=42)
    >>> print(f"Nodes: {dba.N}, Edges: {dba.num_edges}")
    """

    def __init__(
        self,
        n: int,
        m1: int,
        m2: int,
        p: float = 0.5,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ):
        # Validate inputs
        max_m = max(m1, m2)
        if max_m >= n:
            raise ValueError(f"max(m1, m2) ({max_m}) must be < n ({n})")
        if m1 < 1 or m2 < 1:
            raise ValueError(f"m1 and m2 must be >= 1, got m1={m1}, m2={m2}")
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.m1 = m1
        self.m2 = m2
        self.p = p

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph using C++ backend
        G = self._generate_graph(seed or 0)

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed, **kwargs)


    def _generate_graph(self, seed: int) -> gt.Graph:
        """Generate Dual BA graph using C++ extension."""
        G = gt.Graph(directed=False)
        create_dual_barabasi_albert(G, self.n, self.m1, self.m2, self.p, seed)

        # Add sign property (all positive initially)
        # Use vectorized initialization for speed
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G

    def get_expected_degree(self) -> float:
        """Return expected average degree."""
        expected_m = self.p * self.m1 + (1 - self.p) * self.m2
        return 2 * expected_m

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
        return f"N={self.n}_m1={self.m1}_m2={self.m2}_p={self.p}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"DualBarabasiAlbertGT(n={self.n}, m1={self.m1}, m2={self.m2}, p={self.p}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["DualBarabasiAlbertGT"]
