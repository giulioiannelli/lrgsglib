"""
FullyConnectedGT: graph-tool implementation with C++ backend.

Uses native C++ extension for high-performance complete graph generation.
"""
from __future__ import annotations

from typing import Callable, List, Optional, Tuple, Union

import numpy as np
import graph_tool.all as gt

from ....gt_patches.signed_graph_gt import SignedGraphGT
from .cpp import create_complete_graph


class FullyConnectedGT(SignedGraphGT):
    """
    Fully connected (complete) graph using C++ backend.

    Every node is connected to every other node.
    Total edges: N*(N-1)/2

    Parameters
    ----------
    N : int
        Number of nodes in the complete graph.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.
    with_positions : bool, default False
        Whether to compute and store node positions.
    mode_positions : str or callable, default "circular"
        Layout mode for positions ("circular", "spring", or callable).

    Attributes
    ----------
    n : int
        Number of nodes.
    syshape : int
        System shape (equals N for complete graphs).
    syshapePth : str
        Path string representation.

    Examples
    --------
    >>> fc = FullyConnectedGT(N=50, pflip=0.1, seed=42)
    >>> print(f"Nodes: {fc.N}, Edges: {fc.num_edges}")
    Nodes: 50, Edges: 1225

    Notes
    -----
    Complete graphs have N*(N-1)/2 edges for N nodes.
    """

    def __init__(
        self,
        N: int,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        with_positions: bool = False,
        mode_positions: Union[str, Callable] = "circular",
        **kwargs,
    ):
        # Validate inputs
        if N < 1:
            raise ValueError(f"N must be >= 1, got {N}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self._N = N
        self.n = N  # Alias for compatibility
        self.syshape = N
        self.syshapePth = f"N={N}"
        self.with_positions = with_positions
        self.mode_positions = mode_positions

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph using C++ backend
        G = self._generate_graph()

        # Add positions if requested
        if with_positions:
            self._add_positions(G)

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

        # Apply sign flips if requested
        if pflip > 0:
            self.flip_random_fract_edges()

    def _generate_graph(self) -> gt.Graph:
        """Generate complete graph using C++ extension."""
        G = gt.Graph(directed=False)
        create_complete_graph(G, self._N)

        # Add sign property (all positive initially)
        # Use vectorized initialization for speed
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G

    def _add_positions(self, G: gt.Graph) -> None:
        """Add node positions based on mode_positions."""
        pos = G.new_vertex_property("vector<double>")

        if self.mode_positions == "circular":
            # Circular layout
            angles = np.linspace(0, 2 * np.pi, self._N, endpoint=False)
            for i, v in enumerate(G.vertices()):
                pos[v] = [np.cos(angles[i]), np.sin(angles[i])]
        elif self.mode_positions == "spring":
            # Use graph-tool's spring layout
            pos = gt.sfdp_layout(G)
        elif callable(self.mode_positions):
            # Custom layout function
            positions = self.mode_positions(G)
            for v in G.vertices():
                pos[v] = positions[int(v)]
        else:
            raise ValueError(f"Unsupported mode_positions: {self.mode_positions}")

        G.vertex_properties["pos"] = pos

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for this complete graph."""
        return self._N

    def get_expected_num_edges(self) -> int:
        """Return the expected number of edges for this complete graph."""
        return self._N * (self._N - 1) // 2

    def get_degree_distribution(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the degree distribution.

        For a complete graph, all nodes have degree N-1.

        Returns
        -------
        tuple of (degrees, counts)
            Unique degree values and their frequencies.
        """
        # All nodes have degree N-1
        return np.array([self._N - 1]), np.array([self._N])

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"N={self._N}"

    @syshapePth.setter
    def syshapePth(self, value: str) -> None:
        self._syshapePth = value

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"FullyConnectedGT(N={self._N}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["FullyConnectedGT"]
