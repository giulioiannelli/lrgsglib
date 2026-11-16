"""
kRegularGraphGT: graph-tool random k-regular graph.

``method='native'`` uses graph-tool's built-in random_graph with a constant
degree sampler. ``method='pairing'`` uses the exact uniform pairing sampler
shared with the NX engine, so the same seed gives the same graph on both.
The pairing draw takes its own child of the seed (``pairing_stream``) and
leaves the process-global numpy and graph-tool RNGs alone.
"""

from __future__ import annotations

from typing import Optional, Tuple

import graph_tool.all as gt
import graph_tool.generation as gen
import numpy as np

from ....config.const import KREG_METHOD, KREG_METHOD_PAIRING, KREG_METHODS
from ..._shared._regular_pairing import pairing_stream, random_regular_edges
from ..SignedGraphGT import SignedGraphGT


class kRegularGraphGT(SignedGraphGT):
    """
    Random k-regular graph on graph-tool.

    Every node has exactly degree k. The topology comes from graph-tool's
    random_graph with a constant degree sampler (``method='native'``) or
    from the exact uniform pairing sampler (``method='pairing'``).

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Degree of each node. Must satisfy n*k even.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility. ``method='native'`` seeds the
        global numpy and graph-tool RNGs with it; ``method='pairing'``
        draws from ``pairing_stream(seed)`` only.
    method : {'native', 'pairing'}, default 'native'
        Topology sampler: graph-tool's configuration model, or the exact
        uniform pairing sampler (same graph as NX for the same seed).

    Attributes
    ----------
    n : int
        Number of nodes.
    k : int
        Degree of each node.
    syshapePth : str
        Path string, ``N=<n>_k=<k>``, plus ``_method=<method>`` for a
        non-default method.

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
        method: str = KREG_METHOD,
        **kwargs,
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
        if method not in KREG_METHODS:
            raise ValueError(f"method must be one of {KREG_METHODS}")

        self.n = n
        self.k = k
        self.method = method
        self._pairing_seed = pairing_stream(seed)
        # Set before super().__init__ so the lazy path init uses it
        self._syshapePth = f"N={n}_k={k}"
        if method != KREG_METHOD:
            self._syshapePth += f"_method={method}"

        # The native generator draws from the global RNGs
        if seed is not None and method != KREG_METHOD_PAIRING:
            np.random.seed(seed)
            gt.seed_rng(seed)

        G = self._generate_graph()

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed, **kwargs)

    def _generate_graph(self) -> gt.Graph:
        """Generate the k-regular graph with the selected ``method``."""
        if self.method == KREG_METHOD_PAIRING:
            G = self._graph_from_pairing()
        else:
            # graph-tool's random_graph with a constant degree sampler
            # (optimized C++ under the hood)
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

    def _graph_from_pairing(self) -> gt.Graph:
        """Undirected GT graph of the exact pairing sampler's edge list.

        Drawn from ``pairing_stream(seed)``, so the edges equal
        ``random_regular_edges(n, k, pairing_stream(seed))``.
        """
        edges = random_regular_edges(self.n, self.k, self._pairing_seed)
        G = gt.Graph(directed=False)  # noqa: raw-graph
        G.add_vertex(self.n)
        G.add_edge_list(edges)
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

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"kRegularGraphGT(n={self.n}, k={self.k}, method={self.method}, "
            f"edges={self.num_edges}, negative={neg})"
        )


__all__ = ["kRegularGraphGT"]
