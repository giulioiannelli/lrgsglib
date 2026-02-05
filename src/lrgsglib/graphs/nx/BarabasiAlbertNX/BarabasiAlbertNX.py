"""
BarabasiAlbertNX: Scale-free graph via preferential attachment.

This module provides the BarabasiAlbertNX class for generating scale-free
networks using the Barabasi-Albert preferential attachment model.

Example
-------
>>> from lrgsglib.graphs.nx import BarabasiAlbertNX
>>> ba = BarabasiAlbertNX(n=100, m=3, pflip=0.1, seed=42)
>>> ba.flip_random_fract_edges()
>>> ba.N
100
"""

import networkx as nx

from ..RandomGraphNX.RandomGraphNX import RandomGraphNX


# Constants for BarabasiAlbert
BA_PHTABB = "ba"
BA_SGPATH = ""
BA_STDFN = ""


class BarabasiAlbertNX(RandomGraphNX):
    """
    Signed Barabasi-Albert scale-free graph.

    Generates a random graph using the Barabasi-Albert preferential
    attachment model. The resulting graph has a power-law degree distribution
    characteristic of many real-world networks.

    Parameters
    ----------
    n : int
        Number of nodes in the final graph.
    m : int
        Number of edges to attach from a new node to existing nodes.
        Must be less than or equal to n.
    sgpathn : str, default BA_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default BA_STDFN
        Filename suffix used in on-disk exports.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraphNX`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Attributes
    ----------
    m : int
        Number of edges attached per new node.
    gamma : float
        Expected degree exponent (~3 for BA model).

    Notes
    -----
    The Barabasi-Albert model generates graphs with:
    - Power-law degree distribution P(k) ~ k^(-gamma), gamma ~ 3
    - Hub nodes with high connectivity
    - Small-world properties (short average path length)

    The graph is always connected by construction.

    ``pflip`` defines how many edges are marked for sign flips, but weights
    are not set negative until you call ``flip_random_fract_edges`` or
    ``flip_sel_edges``.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import BarabasiAlbertNX
    >>> ba = BarabasiAlbertNX(n=500, m=2, pflip=0.15, seed=42)
    >>> ba.flip_random_fract_edges()
    >>> print(f"Nodes: {ba.N}, Edges: {ba.Ne}")
    Nodes: 500, Edges: 996
    """

    gamma = 3.0  # Expected degree exponent for BA model

    def __init__(
        self,
        n: int,
        m: int,
        sgpathn: str = BA_SGPATH,
        stdFnameSFFX: str = BA_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if m > n:
            raise ValueError(f"m ({m}) must be <= n ({n})")

        self.m = m
        self._init_std_fname(stdFnameSFFX)
        sgpath = BA_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,  # BA graphs are connected by construction
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = BA_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate Barabasi-Albert preferential attachment graph."""
        return nx.barabasi_albert_graph(self.n, self.m)

    def _compute_syshapePth(self) -> str:
        """Compute path string for Barabasi-Albert parameters."""
        return f"N={self.n}_m={self.m}"

    def get_expected_degree(self) -> float:
        """
        Return the expected average degree of the graph.

        For BA model with m edges per new node, the expected average
        degree is approximately 2*m for large n.

        Returns
        -------
        float
            Expected average degree.
        """
        return 2 * self.m

    def get_hub_nodes(self, top_k: int = 10) -> list:
        """
        Return the top-k hub nodes by degree.

        Parameters
        ----------
        top_k : int, default 10
            Number of top nodes to return.

        Returns
        -------
        list
            List of (node, degree) tuples sorted by degree descending.
        """
        degrees = sorted(self.G.degree(), key=lambda x: x[1], reverse=True)
        return degrees[:top_k]
