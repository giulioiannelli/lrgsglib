"""
kRegularGraphGT: graph-tool implementation using native GT random_graph.

Uses graph-tool's built-in random_graph with constant degree sampler
for high-performance k-regular graph generation.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import graph_tool.all as gt
import graph_tool.generation as gen

from ..SignedGraphGT import SignedGraphGT


class kRegularGraphGT(SignedGraphGT):
    """
    Random k-regular graph using graph-tool's native random_graph.

    Every node has exactly degree k. Uses graph-tool's random_graph
    function with a constant degree sampler.

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Degree of each node. Must satisfy n*k even.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    n : int
        Number of nodes.
    k : int
        Degree of each node.
    syshapePth : str
        Path string representation.

    Examples
    --------
    >>> kreg = kRegularGraphGT(n=100, k=4, pflip=0.1, seed=42)
    >>> all(v.out_degree() == 4 for v in kreg.G.vertices())
    True

    Notes
    -----
    k-regular graphs are useful for:
    - Baseline comparisons (zero degree variance)
    - Studying pure structural effects without degree heterogeneity
    - Expander graph constructions
    """

    def __init__(
        self,
        n: int,
        k: int,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ):
        # Validate inputs
        if (n * k) % 2 != 0:
            raise ValueError(f"n*k must be even, got n={n}, k={k}, n*k={n*k}")
        if k >= n:
            raise ValueError(f"k must be less than n, got k={k}, n={n}")
        if k < 0:
            raise ValueError(f"k must be non-negative, got k={k}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.k = k

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph using graph-tool's native random_graph
        G = self._generate_graph()

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

        # Apply sign flips if requested
        if pflip > 0:
            self.flip_random_fract_edges()

    def _generate_graph(self) -> gt.Graph:
        """Generate k-regular graph using graph-tool's random_graph."""
        # Use graph-tool's random_graph with constant degree sampler
        # This is highly optimized C++ under the hood
        G = gen.random_graph(
            self.n,
            lambda: self.k,  # Constant degree sampler
            directed=False,
            parallel_edges=False,
            self_loops=False,
            model="configuration",
        )

        # Add sign property (all positive initially)
        # Use vectorized initialization for speed
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G

    def verify_regularity(self) -> bool:
        """Check that all nodes have degree exactly k.

        Returns
        -------
        bool
            True if all vertices have out-degree k.
        """
        return all(v.out_degree() == self.k for v in self.G.vertices())

    def get_expected_degree(self) -> int:
        """Return the expected (and actual) degree of all nodes."""
        return self.k

    def get_expected_num_edges(self) -> int:
        """Return the expected number of edges."""
        return (self.n * self.k) // 2

    def get_degree_distribution(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the degree distribution.

        For a k-regular graph, all nodes have degree k.

        Returns
        -------
        tuple of (degrees, counts)
            Unique degree values and their frequencies.
        """
        # All nodes have degree k
        return np.array([self.k]), np.array([self.n])

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"N={self.n}_k={self.k}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"kRegularGraphGT(n={self.n}, k={self.k}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["kRegularGraphGT"]
