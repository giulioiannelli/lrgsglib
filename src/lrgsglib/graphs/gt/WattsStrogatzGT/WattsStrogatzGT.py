"""
WattsStrogatzGT: graph-tool implementation of Watts-Strogatz small-world graphs.

Mirrors the API of WattsStrogatzNX for easy switching between backends.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

import graph_tool.all as gt

from ..SignedGraphGT import SignedGraphGT


class WattsStrogatzGT(SignedGraphGT):
    """
    Watts-Strogatz small-world graph using graph-tool backend.

    Creates a small-world network by starting with a ring lattice where
    each node is connected to its k nearest neighbors, then rewiring
    edges with probability p.

    Parameters
    ----------
    n : int
        Number of nodes in the graph.
    k : int
        Each node is initially connected to k nearest neighbors in ring.
        Must be even.
    p : float
        Probability of rewiring each edge (0.0 to 1.0).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    n : int
        Number of nodes.
    k : int
        Initial neighbor connectivity.
    p : float
        Rewiring probability.

    Examples
    --------
    >>> ws = WattsStrogatzGT(n=100, k=4, p=0.3, pflip=0.1, seed=42)
    >>> ws.flip_random_fract_edges()
    >>> print(f"Nodes: {ws.N}, Edges: {ws.num_edges}")

    Notes
    -----
    The Watts-Strogatz model interpolates between regular lattices (p=0)
    and random graphs (p=1), producing small-world networks for
    intermediate values of p.

    Unlike the NetworkX version which uses connected_watts_strogatz_graph,
    this implementation may produce disconnected graphs for high p values.
    """

    def __init__(
        self,
        n: int,
        k: int,
        p: float,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ):
        # Validate inputs
        if k >= n:
            raise ValueError(f"k ({k}) must be < n ({n})")
        if k % 2 != 0:
            raise ValueError(f"k must be even, got {k}")
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.k = k
        self.p = p

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate WS graph
        G = self._generate_graph()

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed, **kwargs)

    def _generate_graph(self) -> gt.Graph:
        """Generate Watts-Strogatz small-world graph.

        Algorithm:
        1. Create ring lattice with k/2 neighbors on each side
        2. For each edge, rewire target with probability p
        """
        G = gt.Graph(directed=False)
        G.add_vertex(self.n)

        # Step 1: Create ring lattice
        # Each node i connects to i+1, i+2, ..., i+k/2 (mod n)
        half_k = self.k // 2
        edges = []
        for i in range(self.n):
            for j in range(1, half_k + 1):
                target = (i + j) % self.n
                edges.append((i, target))

        # Add initial edges
        for src, tgt in edges:
            G.add_edge(src, tgt)

        # Step 2: Rewire edges with probability p
        if self.p > 0:
            edges_to_rewire = []
            for e in G.edges():
                if np.random.random() < self.p:
                    edges_to_rewire.append(e)

            # Get existing edges for duplicate checking
            existing = set()
            for e in G.edges():
                s, t = int(e.source()), int(e.target())
                existing.add((min(s, t), max(s, t)))

            for e in edges_to_rewire:
                src = int(e.source())
                old_tgt = int(e.target())

                # Remove old edge from tracking
                old_key = (min(src, old_tgt), max(src, old_tgt))
                existing.discard(old_key)

                # Find new target (not self, not existing neighbor)
                attempts = 0
                max_attempts = self.n * 2
                while attempts < max_attempts:
                    new_tgt = np.random.randint(0, self.n)
                    if new_tgt != src:
                        new_key = (min(src, new_tgt), max(src, new_tgt))
                        if new_key not in existing:
                            # Rewire: remove old, add new
                            G.remove_edge(e)
                            G.add_edge(src, new_tgt)
                            existing.add(new_key)
                            break
                    attempts += 1

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G

    def get_clustering_coefficient(self) -> float:
        """
        Compute the average clustering coefficient.

        Returns
        -------
        float
            Average clustering coefficient over all nodes.
        """
        from graph_tool.clustering import local_clustering

        cc = local_clustering(self.G)
        return float(np.mean(cc.a))

    def get_average_path_length(self) -> float:
        """
        Compute the average shortest path length.

        Returns
        -------
        float
            Average shortest path length (may be inf if disconnected).
        """
        from graph_tool.topology import shortest_distance

        # Sample-based estimation for large graphs
        if self.N > 1000:
            # Sample 100 source nodes
            sources = np.random.choice(self.N, min(100, self.N), replace=False)
            total = 0
            count = 0
            for s in sources:
                dist = shortest_distance(self.G, source=self.G.vertex(s))
                valid = dist.a[dist.a < self.N]  # Exclude unreachable
                total += valid.sum()
                count += len(valid) - 1  # Exclude self
            return total / count if count > 0 else float('inf')
        else:
            # Full computation for small graphs
            total = 0
            count = 0
            for v in self.G.vertices():
                dist = shortest_distance(self.G, source=v)
                valid = dist.a[dist.a < self.N]
                total += valid.sum()
                count += len(valid) - 1
            return total / count if count > 0 else float('inf')

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"N={self.n}_k={self.k}_p={self.p:.3g}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"WattsStrogatzGT(n={self.n}, k={self.k}, p={self.p}, "
            f"edges={self.num_edges}, negative={neg})"
        )

__all__ = ["WattsStrogatzGT"]
