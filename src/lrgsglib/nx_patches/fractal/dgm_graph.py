"""
DGMgraph: Backward-compatible Dorogovtsev-Goltsev-Mendes graph wrapper.

This module provides the DGMgraph class that wraps the FractalGraph base class
while maintaining full backward compatibility with the existing API.

The DGM graph is a deterministic scale-free network with hierarchical structure.
It has exactly 3*2^n - 2 nodes at iteration n, with degree distribution P(k) ~ k^(-3).

Example
-------
>>> from lrgsglib.nx_patches.fractal import DGMgraph
>>> dgm = DGMgraph(n=5, pflip=0.1, seed=42)
>>> dgm.flip_random_fract_edges()
>>> dgm.N
94
"""

from typing import Any

import networkx as nx

from ...config.const import DGM_ONREP, DGM_PHTABB, DGM_SGPATH, DGM_STDFN
from ..DGMgraph.generators_dgm import dorogovtsev_goltsev_mendes_graph_FastPatch
from ..SignedGraph.SignedGraph import SignedGraph
from .fractal_graph import FractalGraph, FractalNwContainerBase


class DGMgraph(FractalGraph):
    """
    Dorogovtsev-Goltsev-Mendes graph as a SignedGraph.

    The DGM graph is constructed iteratively:
    - Start with a triangle (n=0)
    - At each iteration, add a new node connected to both endpoints of every edge

    This produces a deterministic scale-free network with:
    - Node count: (3^n + 3) / 2
    - Edge count: 3^n
    - Degree exponent: gamma = 3 (exact)
    - Clustering coefficient: non-zero
    - Hierarchical/modular structure

    Parameters
    ----------
    n : int
        Iteration level (number of generations).
    with_positions : bool, default False
        Whether to compute spring layout positions.
    sgpathn : str, default DGM_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default DGM_STDFN
        Filename suffix used in on-disk exports.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraph`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Examples
    --------
    >>> from lrgsglib.nx_patches.fractal import DGMgraph
    >>> dgm = DGMgraph(n=6, pflip=0.1, seed=42)
    >>> dgm.flip_random_fract_edges()
    >>> print(f"Nodes: {dgm.N}, Edges: {dgm.Ne}")
    Nodes: 190, Edges: 1092
    """

    def __init__(
        self,
        n: int,
        *,
        with_positions: bool = False,
        sgpathn: str = DGM_SGPATH,
        stdFnameSFFX: str = DGM_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self._init_std_fname(stdFnameSFFX)
        sgpath = DGM_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            with_positions=with_positions,
            only_const_mode=only_const_mode,
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = DGM_PHTABB + suffix

    def _compute_expected_nodes(self) -> int:
        """Compute expected number of nodes: (3^n + 3) / 2."""
        return (3 ** self.n + 3) // 2

    def _init_network(self) -> None:
        """Build the DGM graph at iteration level n."""
        self.H = dorogovtsev_goltsev_mendes_graph_FastPatch(self.n)
        self.G = nx.convert_node_labels_to_integers(self.H)
        if self.with_positions:
            pos = nx.spring_layout(self.H)
            nx.set_node_attributes(self.H, pos, "pos")
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={self.syshape}"

    def get_central_edge(self, on_g: str = DGM_ONREP) -> tuple:
        """
        Return the central edge of the DGM graph.

        The central edge connects the two original nodes from the
        initial triangle.

        Parameters
        ----------
        on_g : str, default 'G'
            Graph representation to use.

        Returns
        -------
        tuple
            Edge as (node1, node2).
        """
        return (0, 1)

    class nwContainer(FractalNwContainerBase):
        """Network container for DGM graph patterns."""

        def __init__(
            self,
            dgm: SignedGraph,
            iterable: list = None,
            constant: Any = None,
            **kwargs,
        ):
            if iterable is None:
                iterable = []
            # Use 'dgm' as the graph reference for backward compatibility
            self.dgm = dgm
            super().__init__(dgm, iterable, constant, **kwargs)
