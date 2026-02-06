"""
StochasticBlockModelNX: Community-structured random graph.

This module provides the StochasticBlockModelNX class for generating
random graphs with planted community structure using the stochastic
block model (SBM).

Example
-------
>>> from lrgsglib.graphs.nx import StochasticBlockModelNX
>>> sizes = [30, 30, 40]  # Three communities
>>> p_matrix = [[0.3, 0.05, 0.05], [0.05, 0.3, 0.05], [0.05, 0.05, 0.3]]
>>> sbm = StochasticBlockModelNX(sizes=sizes, p_matrix=p_matrix, pflip=0.1)
>>> sbm.flip_random_fract_edges()
"""

from typing import Sequence

import networkx as nx
import numpy as np

from ..RandomGraphNX.RandomGraphNX import RandomGraphNX


# Constants for StochasticBlockModel
SBM_PHTABB = "sbm"
SBM_SGPATH = ""
SBM_STDFN = ""


class StochasticBlockModelNX(RandomGraphNX):
    """
    Signed Stochastic Block Model graph with planted communities.

    Generates a random graph using the stochastic block model where
    nodes are partitioned into communities with different intra- and
    inter-community edge probabilities.

    Parameters
    ----------
    sizes : Sequence[int]
        Sizes of each community (block). Sum determines total nodes.
    p_matrix : Sequence[Sequence[float]]
        Probability matrix where p_matrix[i][j] is the probability
        of an edge between nodes in community i and community j.
        Should be symmetric for undirected graphs.
    sgpathn : str, default SBM_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default SBM_STDFN
        Filename suffix used in on-disk exports.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraphNX`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Attributes
    ----------
    sizes : list[int]
        Community sizes.
    p_matrix : list[list[float]]
        Edge probability matrix.
    num_communities : int
        Number of communities (blocks).
    community_labels : dict
        Mapping from node to community index.

    Notes
    -----
    The stochastic block model is useful for:
    - Benchmarking community detection algorithms
    - Studying spectral properties of modular networks
    - Analyzing dynamics on networks with community structure

    The graph may be disconnected if inter-community probabilities are low.
    By default, the giant component is extracted.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import StochasticBlockModelNX
    >>> # Create graph with 3 equal communities
    >>> sizes = [50, 50, 50]
    >>> p_in, p_out = 0.2, 0.02  # Dense within, sparse between
    >>> p_matrix = [[p_in, p_out, p_out],
    ...             [p_out, p_in, p_out],
    ...             [p_out, p_out, p_in]]
    >>> sbm = StochasticBlockModelNX(sizes=sizes, p_matrix=p_matrix, seed=42)
    >>> print(f"Nodes: {sbm.N}, Communities: {sbm.num_communities}")
    """

    def __init__(
        self,
        sizes: Sequence[int],
        p_matrix: Sequence[Sequence[float]],
        sgpathn: str = SBM_SGPATH,
        stdFnameSFFX: str = SBM_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.sizes = list(sizes)
        self.p_matrix = [list(row) for row in p_matrix]
        self.num_communities = len(sizes)

        # Validate inputs
        self._validate_inputs()

        self._init_std_fname(stdFnameSFFX)
        sgpath = SBM_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=sum(sizes),
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=True,
            **kwargs,
        )

        # Store community labels after graph is built
        if not only_const_mode:
            self._assign_community_labels()

    def _validate_inputs(self) -> None:
        """Validate sizes and p_matrix inputs."""
        if len(self.p_matrix) != self.num_communities:
            raise ValueError(
                f"p_matrix rows ({len(self.p_matrix)}) must match "
                f"number of communities ({self.num_communities})"
            )
        for i, row in enumerate(self.p_matrix):
            if len(row) != self.num_communities:
                raise ValueError(
                    f"p_matrix row {i} has {len(row)} columns, "
                    f"expected {self.num_communities}"
                )
            for j, p in enumerate(row):
                if not 0 <= p <= 1:
                    raise ValueError(
                        f"p_matrix[{i}][{j}] = {p} must be in [0, 1]"
                    )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = SBM_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate stochastic block model graph."""
        return nx.stochastic_block_model(self.sizes, self.p_matrix)

    def _compute_syshapePth(self) -> str:
        """Compute path string for SBM parameters."""
        n = sum(self.sizes)
        k = self.num_communities
        # Compute average intra/inter probabilities for summary
        p_in = np.mean([self.p_matrix[i][i] for i in range(k)])
        p_out_vals = [
            self.p_matrix[i][j]
            for i in range(k)
            for j in range(k)
            if i != j
        ]
        p_out = np.mean(p_out_vals) if p_out_vals else 0
        return f"N={n}_k={k}_pin={p_in:.2g}_pout={p_out:.2g}"

    def _assign_community_labels(self) -> None:
        """Assign community labels to nodes based on original block assignments."""
        self.community_labels = {}
        node_idx = 0
        for comm_idx, size in enumerate(self.sizes):
            for _ in range(size):
                if node_idx < self.G.number_of_nodes():
                    self.community_labels[node_idx] = comm_idx
                node_idx += 1

    def get_community_sizes(self) -> list[int]:
        """
        Return the sizes of each community.

        Returns
        -------
        list[int]
            List of community sizes.
        """
        return self.sizes.copy()

    def get_intra_edge_probability(self) -> float:
        """
        Return the average intra-community edge probability.

        Returns
        -------
        float
            Average probability of edges within communities.
        """
        return float(np.mean([self.p_matrix[i][i] for i in range(self.num_communities)]))

    def get_inter_edge_probability(self) -> float:
        """
        Return the average inter-community edge probability.

        Returns
        -------
        float
            Average probability of edges between communities.
        """
        p_out_vals = [
            self.p_matrix[i][j]
            for i in range(self.num_communities)
            for j in range(self.num_communities)
            if i != j
        ]
        return float(np.mean(p_out_vals)) if p_out_vals else 0.0

    def get_modularity_parameter(self) -> float:
        """
        Return the effective modularity parameter (p_in - p_out) / (p_in + p_out).

        A value close to 1 indicates strong community structure,
        while 0 indicates no structure.

        Returns
        -------
        float
            Modularity parameter in [-1, 1].
        """
        p_in = self.get_intra_edge_probability()
        p_out = self.get_inter_edge_probability()
        total = p_in + p_out
        if total == 0:
            return 0.0
        return (p_in - p_out) / total
