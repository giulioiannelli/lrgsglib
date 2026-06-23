"""
ErdosRenyiGT: graph-tool implementation of Erdos-Renyi random graphs.

Mirrors the API of ErdosRenyiNX for easy switching between backends.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

import graph_tool.all as gt
import graph_tool.generation as gen
import graph_tool.topology as topo

from ....config.const import SG_INIT_NW_DICT
from ..SignedGraphGT import SignedGraphGT
from ..._shared._nw_container import GTnwContainer


class ErdosRenyiGT(SignedGraphGT):
    """
    Erdos-Renyi random graph using graph-tool backend.

    Creates a random graph where each possible edge is included independently
    with probability p. Optionally extracts the giant (largest) connected
    component.

    Parameters
    ----------
    n : int
        Number of nodes.
    p : float
        Probability of edge creation (0.0 to 1.0).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    extract_giant_component : bool
        If True, keep only the largest connected component.
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    n : int
        Original requested number of nodes.
    p : float
        Edge probability.
    N : int
        Actual number of nodes (may be less than n if GC extracted).

    Examples
    --------
    >>> er = ErdosRenyiGT(n=200, p=0.05, pflip=0.2, seed=42)
    >>> er.flip_random_fract_edges()
    >>> print(f"Nodes: {er.N}, Edges: {er.num_edges}")

    Notes
    -----
    With `extract_giant_component=True`, the resulting graph will typically
    have fewer nodes than requested, as only the largest connected component
    is retained.
    """

    # Negative-link (nwDict) pattern container: ErdosRenyi has no central edge,
    # so only the random patterns are built (rand / randXERR), mirroring NX.
    nwContainer = GTnwContainer

    def __init__(
        self,
        n: int,
        p: float,
        pflip: float = 0.0,
        extract_giant_component: bool = True,
        seed: Optional[int] = None,
        init_nw_dict: bool = SG_INIT_NW_DICT,
    ):
        # Validate inputs
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.n = n
        self.p = p
        self.extract_giant_component = extract_giant_component

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Generate Erdos-Renyi graph
        # graph-tool uses a different parametrization
        # We need to compute expected number of edges: E = n*(n-1)/2 * p
        G = self._generate_graph()

        # Extract giant component if requested
        if extract_giant_component:
            G = self._extract_giant_component(G)

        # Set syshapePth before super().__init__ so lazy path init uses it
        self._syshapePth = f"N={G.num_vertices()}_p={p:.3g}"

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed,
                         sgpathn="erdos_renyi_gt",
                         init_nw_dict=init_nw_dict)

    def _generate_graph(self) -> gt.Graph:
        """Generate Erdos-Renyi random graph using graph-tool."""
        G = gt.Graph(directed=False)
        G.add_vertex(self.n)

        # Generate edges: include each possible edge with probability p
        # For efficiency, generate expected number of edges and sample
        max_edges = self.n * (self.n - 1) // 2
        expected_edges = int(max_edges * self.p)

        # For small graphs or high p, iterate all pairs
        # For large sparse graphs, sample edges
        if self.p > 0.5 or self.n < 100:
            # Dense case: iterate all pairs
            for i in range(self.n):
                for j in range(i + 1, self.n):
                    if np.random.random() < self.p:
                        G.add_edge(i, j)
        else:
            # Sparse case: sample approximately correct number of edges
            # Generate random edge indices
            edges_added = set()
            target_edges = expected_edges

            # Add some extra to account for duplicates
            attempts = 0
            max_attempts = target_edges * 10

            while len(edges_added) < target_edges and attempts < max_attempts:
                i = np.random.randint(0, self.n)
                j = np.random.randint(0, self.n)
                if i != j:
                    edge = (min(i, j), max(i, j))
                    if edge not in edges_added:
                        edges_added.add(edge)
                        G.add_edge(edge[0], edge[1])
                attempts += 1

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G

    def _extract_giant_component(self, G: gt.Graph) -> gt.Graph:
        """Extract the largest connected component."""
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

        # Re-add sign property to new graph
        if "sign" not in G_gc.edge_properties:
            sign = G_gc.new_edge_property("int")
            for e in G_gc.edges():
                sign[e] = 1
            G_gc.edge_properties["sign"] = sign

        return G_gc

    @property
    def syshapePth(self) -> str:
        """Parameter string for file paths."""
        return f"N={self.n}_p={self.p:.3g}"

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"ErdosRenyiGT(n={self.n}, p={self.p}, N={self.N}, "
            f"edges={self.num_edges}, negative={neg})"
        )

__all__ = ["ErdosRenyiGT"]
