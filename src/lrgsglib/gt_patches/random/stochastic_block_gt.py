"""
StochasticBlockModelGT: graph-tool implementation of Stochastic Block Model graphs.

Mirrors the API of nx_patches.random.StochasticBlockModel for easy switching between backends.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

import graph_tool.all as gt
import graph_tool.generation as gen
import graph_tool.topology as topo

from ..signed_graph_gt import SignedGraphGT


class StochasticBlockModelGT(SignedGraphGT):
    """
    Stochastic Block Model graph using graph-tool backend.

    Creates a random graph with planted community structure where nodes
    are partitioned into blocks with different intra- and inter-block
    edge probabilities.

    Parameters
    ----------
    sizes : Sequence[int]
        Sizes of each community (block). Sum determines total nodes.
    p_matrix : Sequence[Sequence[float]]
        Probability matrix where p_matrix[i][j] is the probability
        of an edge between nodes in block i and block j.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    extract_giant_component : bool
        If True, keep only the largest connected component.
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    sizes : list[int]
        Community sizes.
    p_matrix : list[list[float]]
        Edge probability matrix.
    num_communities : int
        Number of communities (blocks).
    block_membership : numpy.ndarray
        Block assignment for each node.

    Examples
    --------
    >>> sizes = [30, 30, 40]
    >>> p_matrix = [[0.3, 0.05, 0.05],
    ...             [0.05, 0.3, 0.05],
    ...             [0.05, 0.05, 0.3]]
    >>> sbm = StochasticBlockModelGT(sizes=sizes, p_matrix=p_matrix, seed=42)
    >>> sbm.flip_random_fract_edges()
    >>> print(f"Nodes: {sbm.N}, Communities: {sbm.num_communities}")

    Notes
    -----
    Uses graph-tool's generate_sbm for efficient graph generation.
    The stochastic block model is useful for:
    - Benchmarking community detection algorithms
    - Studying spectral properties of modular networks
    - Analyzing dynamics on networks with community structure
    """

    def __init__(
        self,
        sizes: Sequence[int],
        p_matrix: Sequence[Sequence[float]],
        pflip: float = 0.0,
        extract_giant_component: bool = True,
        seed: Optional[int] = None,
    ):
        self.sizes = list(sizes)
        self.p_matrix = [list(row) for row in p_matrix]
        self.num_communities = len(sizes)
        self.extract_giant_component = extract_giant_component

        # Validate inputs
        self._validate_inputs()

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate SBM graph
        G, self.block_membership = self._generate_graph()

        # Extract giant component if requested
        if extract_giant_component:
            G, self.block_membership = self._extract_giant_component(
                G, self.block_membership
            )

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

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

    def _generate_graph(self) -> tuple[gt.Graph, np.ndarray]:
        """Generate SBM graph using graph-tool.

        Returns
        -------
        tuple[gt.Graph, np.ndarray]
            The graph and block membership array.
        """
        n = sum(self.sizes)
        G = gt.Graph(directed=False)
        G.add_vertex(n)

        # Create block membership array
        block_membership = np.zeros(n, dtype=np.int32)
        idx = 0
        for block_id, size in enumerate(self.sizes):
            block_membership[idx : idx + size] = block_id
            idx += size

        # Generate edges based on probability matrix
        for i in range(n):
            block_i = block_membership[i]
            for j in range(i + 1, n):
                block_j = block_membership[j]
                p = self.p_matrix[block_i][block_j]
                if np.random.random() < p:
                    G.add_edge(i, j)

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G, block_membership

    def _extract_giant_component(
        self, G: gt.Graph, block_membership: np.ndarray
    ) -> tuple[gt.Graph, np.ndarray]:
        """Extract the largest connected component.

        Parameters
        ----------
        G : gt.Graph
            Input graph.
        block_membership : np.ndarray
            Block assignments for original graph.

        Returns
        -------
        tuple[gt.Graph, np.ndarray]
            Giant component graph and updated block membership.
        """
        # Get component labels
        comp, hist = topo.label_components(G)

        # Find largest component
        largest_comp = np.argmax(hist)

        # Create vertex filter for largest component
        vfilt = G.new_vertex_property("bool")
        for v in G.vertices():
            vfilt[v] = comp[v] == largest_comp

        # Create filtered view and copy to new graph
        G.set_vertex_filter(vfilt)
        G_gc = gt.Graph(G, prune=True)
        G.clear_filters()

        # Update block membership for remaining nodes
        old_to_new = {}
        new_idx = 0
        for v in G.vertices():
            if vfilt[v]:
                old_to_new[int(v)] = new_idx
                new_idx += 1

        new_membership = np.array(
            [block_membership[old] for old in sorted(old_to_new.keys())]
        )

        # Re-add sign property to new graph
        if "sign" not in G_gc.edge_properties:
            sign = G_gc.new_edge_property("int")
            for e in G_gc.edges():
                sign[e] = 1
            G_gc.edge_properties["sign"] = sign

        return G_gc, new_membership

    def get_community_sizes(self) -> list[int]:
        """
        Return the actual sizes of each community in the final graph.

        Returns
        -------
        list[int]
            List of community sizes (may differ from input if GC extracted).
        """
        unique, counts = np.unique(self.block_membership, return_counts=True)
        result = [0] * self.num_communities
        for block_id, count in zip(unique, counts):
            if block_id < self.num_communities:
                result[block_id] = count
        return result

    def get_intra_edge_probability(self) -> float:
        """
        Return the average intra-community edge probability.

        Returns
        -------
        float
            Average probability of edges within communities.
        """
        return float(
            np.mean([self.p_matrix[i][i] for i in range(self.num_communities)])
        )

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

    def get_block_membership(self) -> np.ndarray:
        """
        Return the block membership for each node.

        Returns
        -------
        np.ndarray
            Array of block indices for each node.
        """
        return self.block_membership.copy()

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        n = sum(self.sizes)
        k = self.num_communities
        p_in = self.get_intra_edge_probability()
        p_out = self.get_inter_edge_probability()
        return f"N={n}_k={k}_pin={p_in:.2g}_pout={p_out:.2g}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"StochasticBlockModelGT(n={self.N}, communities={self.num_communities}, "
            f"edges={self.num_edges}, negative={neg})"
        )
