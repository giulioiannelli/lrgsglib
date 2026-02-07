"""
RandomGeometricGT: graph-tool implementation using native GT geometric_graph.

Uses graph-tool's built-in geometric_graph for high-performance
spatial graph generation with KD-tree acceleration.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import graph_tool.all as gt
import graph_tool.generation as gen

from ....gt_patches.signed_graph_gt import SignedGraphGT


class RandomGeometricGT(SignedGraphGT):
    """
    Random geometric graph using graph-tool's native geometric_graph.

    Nodes are placed uniformly at random in a unit hypercube [0,1]^dim.
    Two nodes are connected if their Euclidean distance is less than
    the given radius.

    Parameters
    ----------
    n : int
        Number of nodes.
    radius : float
        Connection threshold distance.
    dim : int, default 2
        Dimension of the embedding space.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    n : int
        Number of nodes.
    radius : float
        Connection threshold.
    dim : int
        Embedding dimension.

    Examples
    --------
    >>> rgg = RandomGeometricGT(n=200, radius=0.15, dim=2, pflip=0.1, seed=42)

    Notes
    -----
    Random geometric graphs are useful for:
    - Modeling wireless networks and sensor networks
    - Studying percolation transitions
    - Spatial network analysis
    - Graphs with natural clustering

    The critical radius for connectivity in 2D is approximately
    r_c ~ sqrt(log(n) / (pi * n)).
    """

    def __init__(
        self,
        n: int,
        radius: float,
        dim: int = 2,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ):
        if n <= 0:
            raise ValueError(f"n must be positive, got {n}")
        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius}")
        if dim < 1:
            raise ValueError(f"dim must be >= 1, got {dim}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.radius = radius
        self.dim = dim

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph using graph-tool's native geometric_graph
        G, self._pos = self._generate_graph()

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

        # Apply sign flips if requested
        if pflip > 0:
            self.flip_random_fract_edges()

    def _generate_graph(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Generate random geometric graph using graph-tool."""
        # Generate random positions in unit hypercube
        points = np.random.random((self.n, self.dim))

        # Use graph-tool's geometric_graph (uses KD-tree internally)
        G, pos = gen.geometric_graph(points, self.radius)

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G, pos

    def get_positions(self) -> np.ndarray:
        """
        Return node positions as a numpy array.

        Returns
        -------
        np.ndarray
            Array of shape (n, dim) with positions.
        """
        return np.array([self._pos[v] for v in self.G.vertices()])

    def get_critical_radius(self) -> float:
        """
        Return the approximate critical radius for connectivity in 2D.

        Returns
        -------
        float
            Critical radius r_c ~ sqrt(log(n) / (pi * n)).
        """
        if self.dim != 2:
            raise NotImplementedError("Critical radius formula only for dim=2")
        return float(np.sqrt(np.log(self.n) / (np.pi * self.n)))

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
        return f"N={self.n}_r={self.radius:.3g}_d={self.dim}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"RandomGeometricGT(n={self.n}, radius={self.radius:.3g}, dim={self.dim}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["RandomGeometricGT"]
