"""
Preferential Attachment Variants: Extended scale-free graph models.

This module provides variants of the Barabasi-Albert preferential attachment
model with additional mechanisms like fitness, aging, and initial attractiveness.

Graph Types
-----------
- ExtendedBarabasiAlbert: BA with initial attractiveness parameter
- DualBarabasiAlbert: BA with both preferential attachment and random edges
- HolmeKim: BA with clustering (triad formation)

Example
-------
>>> from lrgsglib.nx_patches.random import ExtendedBarabasiAlbert
>>> g = ExtendedBarabasiAlbert(n=500, m=2, p=0.5, seed=42)  # Add clustering
"""

from typing import Optional

import networkx as nx
import numpy as np

from .random_graph import RandomGraph


# Constants
PA_PHTABB = "pa"
PA_SGPATH = ""
PA_STDFN = ""


class ExtendedBarabasiAlbert(RandomGraph):
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
        Forwarded to SignedGraph (pflip, seed, etc.).

    Notes
    -----
    Attachment probability for node i with degree k_i:
        P(i) = (k_i + a) / sum_j(k_j + a)

    When a=0, this reduces to standard BA.
    As a increases, the degree distribution becomes less heavy-tailed.

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import ExtendedBarabasiAlbert
    >>> # Standard BA-like (a=0 approximation)
    >>> g1 = ExtendedBarabasiAlbert(n=500, m=2, a=0.1, seed=42)
    >>> # More homogeneous (high attractiveness)
    >>> g2 = ExtendedBarabasiAlbert(n=500, m=2, a=10, seed=42)
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
        # Use networkx extended BA graph
        return nx.extended_barabasi_albert_graph(self.n, self.m, p=0, q=0)
        # Note: nx.extended_barabasi_albert_graph uses different extension
        # For true a-parameter, we implement manually below

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

    def _generate_graph(self) -> nx.Graph:
        """Generate extended BA graph."""
        return self._generate_graph_with_attractiveness()

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"N={self.n}_m={self.m}_a={self.a:.2f}"


class DualBarabasiAlbert(RandomGraph):
    """
    Dual Barabasi-Albert model with two attachment modes.

    NetworkX dual BA model where new nodes attach with m1 edges (probability p)
    or m2 edges (probability 1-p) via preferential attachment.

    Parameters
    ----------
    n : int
        Number of nodes.
    m1 : int
        Number of edges in mode 1 (with probability p).
    m2 : int
        Number of edges in mode 2 (with probability 1-p).
    p : float, default 0.5
        Probability of using m1 edges (vs m2 edges).
    sgpathn : str, default PA_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default PA_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraph (pflip, seed, etc.).

    Notes
    -----
    - p=1: All nodes use m1 edges
    - p=0: All nodes use m2 edges
    - 0<p<1: Mixed (heterogeneous attachment)

    This models networks where different nodes have different "budgets"
    for forming connections.

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import DualBarabasiAlbert
    >>> # Half use 2 edges, half use 4 edges
    >>> g = DualBarabasiAlbert(n=500, m1=2, m2=4, p=0.5, seed=42)
    """

    def __init__(
        self,
        n: int,
        m1: int,
        m2: int,
        p: float = 0.5,
        sgpathn: str = PA_SGPATH,
        stdFnameSFFX: str = PA_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if m1 < 1 or m2 < 1:
            raise ValueError(f"m1 and m2 must be >= 1, got m1={m1}, m2={m2}")
        if not 0 <= p <= 1:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if max(m1, m2) >= n:
            raise ValueError(f"max(m1, m2) must be < n")

        self.m1 = m1
        self.m2 = m2
        self.p = p
        self._init_std_fname(stdFnameSFFX)
        sgpath = PA_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = PA_PHTABB + "_dual" + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate dual BA graph."""
        seed = getattr(self, '_rng_seed', None)
        return nx.dual_barabasi_albert_graph(self.n, self.m1, self.m2, self.p, seed=seed)

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"N={self.n}_m1={self.m1}_m2={self.m2}_p={self.p:.2f}"


class HolmeKim(RandomGraph):
    """
    Holme-Kim model: Barabasi-Albert with clustering.

    Extends BA by adding a "triad formation" step. After connecting to a
    node via preferential attachment, with probability p the new node also
    connects to a neighbor of that node, forming triangles.

    Parameters
    ----------
    n : int
        Number of nodes.
    m : int
        Number of edges to attach from each new node.
    p : float, default 0.5
        Probability of triad formation step.
    sgpathn : str, default PA_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default PA_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraph (pflip, seed, etc.).

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

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import HolmeKim
    >>> # High clustering scale-free
    >>> g = HolmeKim(n=500, m=3, p=0.8, seed=42)
    >>> import networkx as nx
    >>> print(f"Clustering: {nx.average_clustering(g.G):.3f}")
    """

    def __init__(
        self,
        n: int,
        m: int,
        p: float = 0.5,
        sgpathn: str = PA_SGPATH,
        stdFnameSFFX: str = PA_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if m < 1:
            raise ValueError(f"m must be >= 1, got {m}")
        if not 0 <= p <= 1:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if m >= n:
            raise ValueError(f"m must be < n, got m={m}, n={n}")

        self.m = m
        self.p = p
        self._init_std_fname(stdFnameSFFX)
        sgpath = PA_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = PA_PHTABB + "_hk" + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate Holme-Kim graph."""
        seed = getattr(self, '_rng_seed', None)
        return nx.powerlaw_cluster_graph(self.n, self.m, self.p, seed=seed)

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"N={self.n}_m={self.m}_p={self.p:.2f}"

    def get_clustering_coefficient(self) -> float:
        """
        Return the average clustering coefficient.

        Returns
        -------
        float
            Average clustering coefficient.
        """
        return nx.average_clustering(self.G)
