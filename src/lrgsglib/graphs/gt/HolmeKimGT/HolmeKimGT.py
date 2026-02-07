"""
HolmeKimGT: graph-tool implementation with C++ backend.

Uses native C++ extension for high-performance Holme-Kim
powerlaw cluster graph generation.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
import graph_tool.all as gt

from ....gt_patches.signed_graph_gt import SignedGraphGT
from .cpp import create_holme_kim


class HolmeKimGT(SignedGraphGT):
    """
    Holme-Kim powerlaw cluster graph using C++ backend.

    Extends Barabasi-Albert with a "triad formation" step. After connecting
    via preferential attachment, with probability p, the new node connects
    to a neighbor of the last target, forming triangles.

    Parameters
    ----------
    n : int
        Number of nodes in final graph.
    m : int
        Number of edges each new node creates.
    p : float
        Probability of triad formation (0=pure BA, 1=max clustering).
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
    p : float
        Triad formation probability.

    Examples
    --------
    >>> hk = HolmeKimGT(n=500, m=3, p=0.8, pflip=0.1, seed=42)
    >>> print(f"Nodes: {hk.N}, Edges: {hk.num_edges}")

    Notes
    -----
    - p=0: Standard BA (no extra clustering)
    - p=1: Maximum clustering (always try to form triangles)

    This model produces scale-free networks with higher clustering
    coefficient than standard BA, more realistic for social networks.

    References
    ----------
    Holme, P., & Kim, B. J. (2002). Growing scale-free networks with
    tunable clustering. Physical Review E, 65(2), 026107.
    """

    def __init__(
        self,
        n: int,
        m: int,
        p: float = 0.5,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ):
        # Validate inputs
        if m >= n:
            raise ValueError(f"m must be < n, got m={m}, n={n}")
        if m < 1:
            raise ValueError(f"m must be >= 1, got {m}")
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.m = m
        self.p = p

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
        """Generate Holme-Kim graph using C++ extension."""
        G = gt.Graph(directed=False)
        create_holme_kim(G, self.n, self.m, self.p, seed)

        # Add sign property (all positive initially)
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

    def get_clustering_coefficient(self) -> float:
        """Return the average local clustering coefficient."""
        return float(gt.global_clustering(self.G)[0])

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"N={self.n}_m={self.m}_p={self.p:.2f}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"HolmeKimGT(n={self.n}, m={self.m}, p={self.p:.2f}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["HolmeKimGT"]
