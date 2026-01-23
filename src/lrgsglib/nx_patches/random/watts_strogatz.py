"""
WattsStrogatz: Backward-compatible Watts-Strogatz small-world graph wrapper.

This module provides the WattsStrogatz class that wraps the RandomGraph base class
while maintaining full backward compatibility with the existing API.

Example
-------
>>> from lrgsglib.nx_patches.random import WattsStrogatz
>>> ws = WattsStrogatz(n=100, k=4, p=0.3, pflip=0.1, seed=42)
>>> ws.flip_random_fract_edges()
>>> ws.N
100
"""

import networkx as nx

from ...config.const import WS_PHTABB, WS_SGPATH, WS_STDFN
from .random_graph import RandomGraph


class WattsStrogatz(RandomGraph):
    """
    Signed Watts-Strogatz small-world graph.

    Parameters
    ----------
    n : int
        Number of nodes in the graph.
    k : int
        Each node is joined with its k nearest neighbors in a ring topology.
    p : float
        Probability of rewiring each edge.
    sgpathn : str, default WS_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default WS_STDFN
        Filename suffix used in on-disk exports.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraph`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Notes
    -----
    Uses ``nx.connected_watts_strogatz_graph`` to ensure the resulting
    graph is connected.

    ``pflip`` defines how many edges are marked for sign flips, but weights
    are not set negative until you call ``flip_random_fract_edges`` or
    ``flip_sel_edges``, unless weights already exist on the input graph.

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import WattsStrogatz
    >>> ws = WattsStrogatz(n=100, k=4, p=0.3, pflip=0.1, seed=1)
    >>> ws.flip_random_fract_edges()
    """

    def __init__(
        self,
        n: int,
        k: int,
        p: float,
        sgpathn: str = WS_SGPATH,
        stdFnameSFFX: str = WS_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.k = k
        self.p = p
        self._init_std_fname(stdFnameSFFX)
        sgpath = WS_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,  # connected_watts_strogatz is already connected
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = WS_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate connected Watts-Strogatz small-world graph."""
        return nx.connected_watts_strogatz_graph(self.n, self.k, self.p)

    def _compute_syshapePth(self) -> str:
        """Compute path string for Watts-Strogatz parameters."""
        return f"N={self.n}_p={self.p:.3g}"
