"""
BipartiteGraphGT: graph-tool implementation of bipartite (two-mode) networks.

Uses native graph-tool operations for high-performance bipartite graph generation.
"""
from __future__ import annotations

from typing import Optional, Set, Tuple

import numpy as np
import graph_tool.all as gt

from ....gt_patches.signed_graph_gt import SignedGraphGT


class BipartiteGraphGT(SignedGraphGT):
    """
    Random bipartite (two-mode) graph using graph-tool.

    Creates a bipartite graph with two disjoint node sets where edges
    only connect nodes from different sets. Edge probability is uniform.

    Parameters
    ----------
    n1 : int
        Number of nodes in the first partition (top nodes).
    n2 : int
        Number of nodes in the second partition (bottom nodes).
    p : float, default 0.1
        Probability of edge between any top-bottom node pair.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    n1 : int
        Size of first partition.
    n2 : int
        Size of second partition.
    p : float
        Edge probability.
    partition : VertexPropertyMap
        Vertex property marking partition (0 = top, 1 = bottom).

    Examples
    --------
    >>> bg = BipartiteGraphGT(n1=50, n2=100, p=0.05, seed=42)
    >>> print(f"Top: {bg.n1}, Bottom: {bg.n2}, Edges: {bg.num_edges}")

    Notes
    -----
    Bipartite graphs have special spectral properties:
    - Eigenvalues are symmetric around zero
    - The spectrum encodes the singular values of the biadjacency matrix
    - Useful for studying two-mode networks and recommendation systems
    """

    def __init__(
        self,
        n1: int,
        n2: int,
        p: float = 0.1,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ):
        if n1 < 1 or n2 < 1:
            raise ValueError(f"n1 and n2 must be >= 1, got n1={n1}, n2={n2}")
        if not 0 <= p <= 1:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n1 = n1
        self.n2 = n2
        self.p = p
        self.syshape = (n1, n2)

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph
        G, self.partition = self._generate_graph()

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

    def _generate_graph(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Generate random bipartite graph."""
        G = gt.Graph(directed=False)
        n_total = self.n1 + self.n2

        # Add vertices
        G.add_vertex(n_total)

        # Create partition property
        partition = G.new_vertex_property("int")
        for v in range(self.n1):
            partition[v] = 0  # top nodes
        for v in range(self.n1, n_total):
            partition[v] = 1  # bottom nodes
        G.vertex_properties["bipartite"] = partition

        # Generate edges using vectorized random sampling
        # For large graphs, sample expected number of edges directly
        expected_edges = int(self.n1 * self.n2 * self.p)

        if self.p < 0.5 and expected_edges < self.n1 * self.n2 * 0.3:
            # Sparse case: sample edges directly (faster)
            num_possible = self.n1 * self.n2
            num_edges = np.random.binomial(num_possible, self.p)

            if num_edges > 0:
                # Sample edge indices without replacement
                edge_indices = np.random.choice(num_possible, size=num_edges, replace=False)

                # Convert to (top, bottom) pairs
                top_nodes = edge_indices // self.n2
                bottom_nodes = edge_indices % self.n2 + self.n1

                # Add edges as batch
                edge_list = np.column_stack((top_nodes, bottom_nodes))
                G.add_edge_list(edge_list)
        else:
            # Dense case: generate all potential edges and filter
            top_idx, bottom_idx = np.meshgrid(
                np.arange(self.n1),
                np.arange(self.n1, self.n1 + self.n2),
                indexing='ij'
            )
            mask = np.random.random((self.n1, self.n2)) < self.p
            edges = np.column_stack((top_idx[mask].ravel(), bottom_idx[mask].ravel()))
            if len(edges) > 0:
                G.add_edge_list(edges)

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G, partition

    def get_top_nodes(self) -> Set[int]:
        """Return set of top partition node indices."""
        return set(range(self.n1))

    def get_bottom_nodes(self) -> Set[int]:
        """Return set of bottom partition node indices."""
        return set(range(self.n1, self.n1 + self.n2))

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes."""
        return self.n1 + self.n2

    def get_expected_num_edges(self) -> int:
        """Return the expected number of edges."""
        return int(self.n1 * self.n2 * self.p)

    def get_biadjacency_matrix(self, sparse: bool = False) -> np.ndarray:
        """
        Return the biadjacency matrix of the bipartite graph.

        The biadjacency matrix B has shape (n1, n2) where B[i,j] = 1
        if there is an edge between top node i and bottom node j.

        Parameters
        ----------
        sparse : bool, default False
            If True, return a scipy sparse matrix.

        Returns
        -------
        ndarray or sparse matrix
            Biadjacency matrix of shape (n1, n2).
        """
        if sparse:
            from scipy import sparse as sp
            rows = []
            cols = []
            for e in self.G.edges():
                s, t = int(e.source()), int(e.target())
                if s < self.n1:
                    rows.append(s)
                    cols.append(t - self.n1)
                else:
                    rows.append(t)
                    cols.append(s - self.n1)
            return sp.csr_matrix(
                (np.ones(len(rows), dtype=int), (rows, cols)),
                shape=(self.n1, self.n2)
            )
        else:
            B = np.zeros((self.n1, self.n2), dtype=int)
            for e in self.G.edges():
                s, t = int(e.source()), int(e.target())
                if s < self.n1:
                    B[s, t - self.n1] = 1
                else:
                    B[t, s - self.n1] = 1
            return B

    def get_density(self) -> float:
        """Return the edge density of the bipartite graph."""
        return self.num_edges / (self.n1 * self.n2)

    def is_bipartite(self) -> bool:
        """Verify that the graph is bipartite."""
        # By construction, always bipartite
        return True

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"n1={self.n1}_n2={self.n2}_p={self.p:.3g}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"BipartiteGraphGT(n1={self.n1}, n2={self.n2}, p={self.p:.3g}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["BipartiteGraphGT"]
