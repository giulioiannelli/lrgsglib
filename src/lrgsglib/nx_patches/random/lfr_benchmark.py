"""
LFRBenchmark: Lancichinetti-Fortunato-Radicchi benchmark graph.

The LFR benchmark generates graphs with known community structure and
power-law degree distribution, commonly used for evaluating community
detection algorithms.

Example
-------
>>> from lrgsglib.nx_patches.random import LFRBenchmark
>>> lfr = LFRBenchmark(n=250, tau1=3, tau2=1.5, mu=0.1, average_degree=5, seed=42)
>>> communities = lfr.get_communities()
>>> print(f"Nodes: {lfr.N}, Communities: {len(communities)}")
"""

from typing import Optional, Sequence

import networkx as nx
import numpy as np

from .random_graph import RandomGraph


# Constants
LFR_PHTABB = "lfr"
LFR_SGPATH = ""
LFR_STDFN = ""


class LFRBenchmark(RandomGraph):
    """
    Signed LFR benchmark graph for community detection.

    The LFR (Lancichinetti-Fortunato-Radicchi) benchmark generates graphs
    with planted community structure and heterogeneous degree/community size
    distributions following power laws.

    Parameters
    ----------
    n : int
        Number of nodes.
    tau1 : float
        Power-law exponent for degree distribution. Typically 2 < tau1 < 3.
    tau2 : float
        Power-law exponent for community size distribution. Typically 1 < tau2 < 2.
    mu : float
        Mixing parameter (fraction of inter-community edges). 0 < mu < 1.
        Lower mu = stronger community structure.
    average_degree : float, optional
        Average node degree. If None, defaults to ~20.
    min_degree : int, optional
        Minimum node degree. If None, computed from average_degree.
    max_degree : int, optional
        Maximum node degree. If None, defaults to n.
    min_community : int, optional
        Minimum community size. Defaults to ~10.
    max_community : int, optional
        Maximum community size. Defaults to ~50.
    sgpathn : str, default LFR_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default LFR_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraph (pflip, seed, etc.).

    Attributes
    ----------
    tau1 : float
        Degree distribution exponent.
    tau2 : float
        Community size distribution exponent.
    mu : float
        Mixing parameter.
    communities : list[set]
        List of node communities.

    Notes
    -----
    The LFR benchmark is useful for:
    - Evaluating community detection algorithms
    - Testing modularity-based methods
    - Studying the detectability threshold

    Key properties:
    - Degree distribution: P(k) ~ k^(-tau1)
    - Community size distribution: P(s) ~ s^(-tau2)
    - Mixing parameter mu controls difficulty (higher = harder to detect)

    References
    ----------
    Lancichinetti, A., Fortunato, S., & Radicchi, F. (2008).
    Benchmark graphs for testing community detection algorithms.
    Physical Review E, 78(4), 046110.

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import LFRBenchmark
    >>> # Create LFR graph with clear communities
    >>> lfr = LFRBenchmark(n=250, tau1=3, tau2=1.5, mu=0.1, average_degree=5, seed=42)
    >>> print(f"Nodes: {lfr.N}, Communities: {len(lfr.get_communities())}")

    >>> # Harder case (more mixing)
    >>> lfr_hard = LFRBenchmark(n=250, tau1=3, tau2=1.5, mu=0.4, average_degree=5, seed=42)
    """

    def __init__(
        self,
        n: int,
        tau1: float,
        tau2: float,
        mu: float,
        average_degree: Optional[float] = None,
        min_degree: Optional[int] = None,
        max_degree: Optional[int] = None,
        min_community: Optional[int] = None,
        max_community: Optional[int] = None,
        sgpathn: str = LFR_SGPATH,
        stdFnameSFFX: str = LFR_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        # Validate parameters
        if tau1 <= 1:
            raise ValueError(f"tau1 must be > 1, got {tau1}")
        if tau2 <= 1:
            raise ValueError(f"tau2 must be > 1, got {tau2}")
        if not 0 < mu < 1:
            raise ValueError(f"mu must be in (0, 1), got {mu}")

        self.tau1 = tau1
        self.tau2 = tau2
        self.mu = mu

        # NetworkX requires exactly one of min_degree or average_degree
        # Prefer average_degree if provided, otherwise use min_degree
        if average_degree is not None and min_degree is not None:
            raise ValueError("Provide only one of average_degree or min_degree, not both")

        self.average_degree = average_degree
        self.min_degree = min_degree

        # If neither provided, default to average_degree (conservative value)
        if self.average_degree is None and self.min_degree is None:
            self.average_degree = 5.0

        self.max_degree = max_degree if max_degree else n

        # Set defaults for community size parameters (LFR is sensitive to these)
        self.min_community = min_community if min_community else 20
        self.max_community = max_community if max_community else max(50, n // 3)

        # Ensure max_community >= min_community
        if self.max_community < self.min_community:
            self.max_community = self.min_community

        self._init_std_fname(stdFnameSFFX)
        sgpath = LFR_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,  # LFR should be connected
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = LFR_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate LFR benchmark graph."""
        # Get seed from kwargs if available
        seed = getattr(self, '_rng_seed', None)

        # Build kwargs - only pass one of average_degree or min_degree
        kwargs = {
            'n': self.n,
            'tau1': self.tau1,
            'tau2': self.tau2,
            'mu': self.mu,
            'max_degree': self.max_degree,
            'min_community': self.min_community,
            'max_community': self.max_community,
            'seed': seed,
        }

        if self.average_degree is not None:
            kwargs['average_degree'] = self.average_degree
        else:
            kwargs['min_degree'] = self.min_degree

        G = nx.generators.community.LFR_benchmark_graph(**kwargs)
        return G

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"N={self.n}_t1={self.tau1:.2f}_t2={self.tau2:.2f}_mu={self.mu:.2f}"

    def get_communities(self) -> list:
        """
        Return the planted communities.

        Returns
        -------
        list[set]
            List of sets, each containing node IDs in a community.
        """
        communities = {frozenset(self.G.nodes[v]["community"]) for v in self.G}
        return [set(c) for c in communities]

    def get_community_sizes(self) -> list:
        """
        Return sizes of planted communities.

        Returns
        -------
        list[int]
            Sorted list of community sizes.
        """
        return sorted([len(c) for c in self.get_communities()], reverse=True)

    def get_modularity(self) -> float:
        """
        Compute modularity of the planted partition.

        Returns
        -------
        float
            Modularity Q of the planted communities.
        """
        communities = self.get_communities()
        return nx.algorithms.community.modularity(self.G, communities)

    def get_node_community(self, node: int) -> set:
        """
        Return the community containing the given node.

        Parameters
        ----------
        node : int
            Node ID.

        Returns
        -------
        set
            Set of nodes in the same community.
        """
        return set(self.G.nodes[node]["community"])

    def get_intra_inter_edges(self) -> tuple:
        """
        Count intra-community and inter-community edges.

        Returns
        -------
        tuple[int, int]
            (intra_edges, inter_edges) counts.
        """
        intra = 0
        inter = 0
        for u, v in self.G.edges():
            comm_u = self.G.nodes[u]["community"]
            comm_v = self.G.nodes[v]["community"]
            if comm_u == comm_v:
                intra += 1
            else:
                inter += 1
        return intra, inter

    def get_actual_mu(self) -> float:
        """
        Compute the actual mixing parameter from the graph.

        Returns
        -------
        float
            Fraction of edges that are inter-community.
        """
        intra, inter = self.get_intra_inter_edges()
        total = intra + inter
        return inter / total if total > 0 else 0.0
