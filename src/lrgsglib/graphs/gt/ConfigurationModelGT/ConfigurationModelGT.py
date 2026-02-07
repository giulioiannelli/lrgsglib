"""
ConfigurationModelGT: graph-tool implementation using native GT random_graph.

Uses graph-tool's built-in random_graph with a degree sequence sampler
for high-performance configuration model graph generation.
"""
from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
import graph_tool.all as gt
import graph_tool.generation as gen

from ....gt_patches.signed_graph_gt import SignedGraphGT


class ConfigurationModelGT(SignedGraphGT):
    """
    Random graph with prescribed degree sequence using graph-tool.

    Uses the configuration model to generate a simple graph (no self-loops
    or multi-edges) with the specified degree sequence.

    Parameters
    ----------
    degree_sequence : Sequence[int]
        Degree sequence. Sum must be even. Length determines n.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    degree_sequence : list[int]
        The prescribed degree sequence.
    n : int
        Number of nodes.
    mean_degree : float
        Mean of the degree sequence.

    Examples
    --------
    >>> degrees = [3, 3, 3, 3, 4, 4, 4, 4, 5, 5]  # sum must be even
    >>> cm = ConfigurationModelGT(degree_sequence=degrees, pflip=0.1, seed=42)
    >>> cm.flip_random_fract_edges()

    Notes
    -----
    The configuration model is useful for:
    - Studying effects of degree heterogeneity
    - Creating null models with matched degree sequences
    - Testing community detection under degree constraints
    """

    def __init__(
        self,
        degree_sequence: Sequence[int],
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ):
        self.degree_sequence = list(degree_sequence)
        self.n = len(self.degree_sequence)

        # Validation
        if sum(self.degree_sequence) % 2 != 0:
            raise ValueError("Sum of degree sequence must be even")
        if any(d < 0 for d in self.degree_sequence):
            raise ValueError("Degree sequence must have non-negative values")
        if any(d >= self.n for d in self.degree_sequence):
            raise ValueError("All degrees must be less than n")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.mean_degree = float(np.mean(self.degree_sequence))

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate graph using graph-tool's native random_graph
        G = self._generate_graph()

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

    def _generate_graph(self) -> gt.Graph:
        """Generate configuration model graph using graph-tool's random_graph."""
        # Create a degree sampler that cycles through the sequence
        deg_iter = iter(self.degree_sequence)

        def deg_sampler():
            try:
                return next(deg_iter)
            except StopIteration:
                # Should never happen if n is correct
                return 0

        # Use graph-tool's random_graph with the degree sampler
        G = gen.random_graph(
            self.n,
            deg_sampler,
            directed=False,
            parallel_edges=False,
            self_loops=False,
            model="configuration",
        )

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        return G

    def get_expected_num_edges(self) -> int:
        """Return the expected number of edges."""
        return sum(self.degree_sequence) // 2

    def get_degree_statistics(self) -> dict:
        """
        Return statistics of the prescribed degree sequence.

        Returns
        -------
        dict
            Dictionary with mean, std, min, max of degree sequence.
        """
        return {
            'mean': float(np.mean(self.degree_sequence)),
            'std': float(np.std(self.degree_sequence)),
            'min': int(np.min(self.degree_sequence)),
            'max': int(np.max(self.degree_sequence)),
        }

    def get_degree_distribution(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the actual degree distribution of the generated graph.

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
        return f"N={self.n}_kavg={self.mean_degree:.2f}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"ConfigurationModelGT(n={self.n}, mean_k={self.mean_degree:.2f}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["ConfigurationModelGT"]
