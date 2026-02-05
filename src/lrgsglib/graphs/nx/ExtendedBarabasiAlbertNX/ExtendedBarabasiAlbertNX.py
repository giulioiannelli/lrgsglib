"""
ExtendedBarabasiAlbertNX: Extended BA model with initial attractiveness.

In the standard BA model, attachment probability is proportional to degree.
This extension adds an initial attractiveness parameter 'a' so that
probability ~ degree + a, allowing new nodes to attract connections.

Example
-------
>>> from lrgsglib.graphs.nx import ExtendedBarabasiAlbertNX
>>> g = ExtendedBarabasiAlbertNX(n=500, m=2, a=1.0, seed=42)
"""

from typing import Optional

import networkx as nx
import numpy as np

from ..RandomGraphNX.RandomGraphNX import RandomGraphNX


# Constants
PA_PHTABB = "pa"
PA_SGPATH = ""
PA_STDFN = ""


class ExtendedBarabasiAlbertNX(RandomGraphNX):
    """
    Extended Barabasi-Albert model with initial attractiveness.

    In the standard BA model, attachment probability is proportional to degree.
    This extension adds an initial attractiveness parameter 'a' so that
    probability ~ degree + a, allowing new nodes to attract connections.

    Parameters
    ----------
    n : int
        Number of nodes.
    m : int
        Number of edges to attach from a new node to existing nodes.
    a : float, default 1.0
        Initial attractiveness. Higher values reduce the "rich get richer"
        effect and lead to more homogeneous degree distributions.
    sgpathn : str, default PA_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default PA_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraphNX (pflip, seed, etc.).

    Notes
    -----
    Attachment probability for node i with degree k_i:
        P(i) = (k_i + a) / sum_j(k_j + a)

    When a=0, this reduces to standard BA.
    As a increases, the degree distribution becomes less heavy-tailed.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import ExtendedBarabasiAlbertNX
    >>> # Standard BA-like (a=0 approximation)
    >>> g1 = ExtendedBarabasiAlbertNX(n=500, m=2, a=0.1, seed=42)
    >>> # More homogeneous (high attractiveness)
    >>> g2 = ExtendedBarabasiAlbertNX(n=500, m=2, a=10, seed=42)
    """

    def __init__(
        self,
        n: int,
        m: int,
        a: float = 1.0,
        sgpathn: str = PA_SGPATH,
        stdFnameSFFX: str = PA_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if m < 1:
            raise ValueError(f"m must be >= 1, got {m}")
        if a < 0:
            raise ValueError(f"a must be >= 0, got {a}")
        if m >= n:
            raise ValueError(f"m must be < n, got m={m}, n={n}")

        self.m = m
        self.a = a
        self._init_std_fname(stdFnameSFFX)
        sgpath = PA_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,  # BA is always connected
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = PA_PHTABB + "_ext" + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate extended BA graph with initial attractiveness."""
        return self._generate_graph_with_attractiveness()

    def _generate_graph_with_attractiveness(self) -> nx.Graph:
        """Generate BA graph with initial attractiveness parameter a."""
        G = nx.Graph()

        # Start with a small complete graph
        initial_nodes = self.m + 1
        G.add_nodes_from(range(initial_nodes))
        for i in range(initial_nodes):
            for j in range(i + 1, initial_nodes):
                G.add_edge(i, j)

        # Get seed for reproducibility
        seed = getattr(self, '_rng_seed', None)
        rng = np.random.default_rng(seed)

        # Add remaining nodes
        for new_node in range(initial_nodes, self.n):
            G.add_node(new_node)

            # Compute attachment probabilities
            degrees = np.array([G.degree(node) for node in range(new_node)])
            probs = degrees + self.a
            probs = probs / probs.sum()

            # Select m nodes to attach to
            targets = rng.choice(
                new_node, size=self.m, replace=False, p=probs
            )

            for target in targets:
                G.add_edge(new_node, target)

        return G

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"N={self.n}_m={self.m}_a={self.a:.2f}"
