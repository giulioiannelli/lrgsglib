"""
BipartiteGraphNX: Signed bipartite (two-mode) network.

Unified class supporting both ER-style random bipartite graphs and
bipartite configuration model graphs from prescribed degree sequences.

Example
-------
>>> from lrgsglib.graphs.nx import BipartiteGraphNX
>>> # Mode 1: ER random bipartite
>>> bg = BipartiteGraphNX(n1=20, n2=30, p=0.1, pflip=0.1, seed=42)
>>> bg.flip_random_fract_edges()
>>> print(f"Nodes: {bg.N}, Edges: {bg.G.number_of_edges()}")
>>>
>>> # Mode 2: Configuration model from degree sequences
>>> bg2 = BipartiteGraphNX(
...     top_degrees=[3, 3, 2, 2],
...     bottom_degrees=[2, 2, 2, 2, 2],
...     seed=42,
... )
>>> print(f"Nodes: {bg2.N}, Mode: {bg2.mode}")
"""

from typing import Optional, Sequence

import networkx as nx
import numpy as np
from networkx.algorithms import bipartite

from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


# Constants
BP_PHTABB = "bp"
BP_SGPATH = ""
BP_STDFN = ""


class BipartiteGraphNX(SignedGraphNX):
    """
    Signed bipartite (two-mode) graph with two build modes.

    Supports two build modes, auto-detected from constructor arguments:

    **Random mode** (``mode="random"``): Creates a bipartite graph with
    edges added independently with probability *p* between partitions.
    Requires ``n1``, ``n2``, and ``p``.

    **Degree-sequence mode** (``mode="degree_sequence"``): Creates a
    bipartite configuration model where top and bottom nodes have
    prescribed degree sequences. Requires ``top_degrees`` and
    ``bottom_degrees``.

    Parameters
    ----------
    n1 : int, optional
        Number of nodes in the first partition (top nodes).
        Required for random mode; inferred from ``top_degrees`` in
        degree-sequence mode.
    n2 : int, optional
        Number of nodes in the second partition (bottom nodes).
        Required for random mode; inferred from ``bottom_degrees`` in
        degree-sequence mode.
    p : float, optional
        Probability of edge between any top-bottom node pair.
        Only used in random mode. Default 0.1.
    top_degrees : Sequence[int], optional
        Degree sequence for top nodes. Triggers degree-sequence mode.
    bottom_degrees : Sequence[int], optional
        Degree sequence for bottom nodes. Triggers degree-sequence mode.
    sgpathn : str, default BP_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default BP_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraphNX (pflip, seed, etc.).

    Attributes
    ----------
    mode : str
        Build mode: ``"random"`` or ``"degree_sequence"``.
    n1 : int
        Size of first partition.
    n2 : int
        Size of second partition.
    p : float
        Edge probability (random mode) or ``0.0`` (degree-sequence mode).
    top_nodes : set
        Nodes in the first partition.
    bottom_nodes : set
        Nodes in the second partition.
    top_degrees : list[int] or None
        Prescribed top degree sequence (degree-sequence mode only).
    bottom_degrees : list[int] or None
        Prescribed bottom degree sequence (degree-sequence mode only).

    Notes
    -----
    Bipartite graphs have special spectral properties:
    - Eigenvalues are symmetric around zero
    - The spectrum encodes the singular values of the biadjacency matrix
    - Useful for studying two-mode networks and recommendation systems

    The bipartite structure is stored as node attributes:
    - ``G.nodes[v]['bipartite']`` is 0 for top nodes, 1 for bottom nodes

    Examples
    --------
    >>> from lrgsglib.graphs.nx import BipartiteGraphNX
    >>> # Random bipartite graph
    >>> bg = BipartiteGraphNX(n1=50, n2=100, p=0.05, seed=42)
    >>> print(f"Top: {len(bg.top_nodes)}, Bottom: {len(bg.bottom_nodes)}")
    Top: 50, Bottom: 100

    >>> # Degree-sequence bipartite graph
    >>> bg2 = BipartiteGraphNX(
    ...     top_degrees=[3, 3, 2, 2],
    ...     bottom_degrees=[2, 2, 2, 2, 2],
    ...     seed=42,
    ... )
    >>> print(f"Mode: {bg2.mode}, Nodes: {bg2.N}")
    """

    def __init__(
        self,
        n1: Optional[int] = None,
        n2: Optional[int] = None,
        p: float = 0.1,
        top_degrees: Optional[Sequence[int]] = None,
        bottom_degrees: Optional[Sequence[int]] = None,
        sgpathn: str = BP_SGPATH,
        stdFnameSFFX: str = BP_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        # --- Resolve build mode ---
        has_deg = (
            top_degrees is not None and bottom_degrees is not None
        )
        has_partial_deg = (
            (top_degrees is not None) != (bottom_degrees is not None)
        )

        if has_partial_deg:
            raise ValueError(
                "Both top_degrees and bottom_degrees must be provided "
                "together for degree-sequence mode."
            )

        if has_deg:
            self.mode: str = "degree_sequence"
            self.top_degrees: Optional[list[int]] = list(top_degrees)
            self.bottom_degrees: Optional[list[int]] = list(
                bottom_degrees
            )

            if sum(self.top_degrees) != sum(self.bottom_degrees):
                raise ValueError(
                    f"Sum of degree sequences must be equal: "
                    f"sum(top)={sum(self.top_degrees)} != "
                    f"sum(bottom)={sum(self.bottom_degrees)}"
                )

            self.n1 = len(self.top_degrees)
            self.n2 = len(self.bottom_degrees)
            self.p = 0.0
        else:
            if n1 is None or n2 is None:
                raise ValueError(
                    "Either (n1, n2) for random mode or "
                    "(top_degrees, bottom_degrees) for "
                    "degree-sequence mode must be provided."
                )
            self.mode = "random"
            self.top_degrees = None
            self.bottom_degrees = None
            self.n1 = n1
            self.n2 = n2
            self.p = p

            if self.n1 < 1 or self.n2 < 1:
                raise ValueError(
                    f"n1 and n2 must be >= 1, got n1={n1}, n2={n2}"
                )
            if not 0 <= self.p <= 1:
                raise ValueError(
                    f"p must be in [0, 1], got {self.p}"
                )

        self.only_const_mode = only_const_mode
        self._init_std_fname(stdFnameSFFX)
        self.sgpathn = BP_PHTABB if not sgpathn else sgpathn
        self.syshape = (self.n1, self.n2)
        self.syshapePth = self._compute_syshapePth()

        if not only_const_mode:
            self._init_bipartite()
        else:
            self.G = nx.Graph()
            self.top_nodes = set()
            self.bottom_nodes = set()

        super().__init__(self.G, **kwargs)

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        if self.mode == "degree_sequence":
            self.std_fname = BP_PHTABB + "_deg" + suffix
        else:
            self.std_fname = BP_PHTABB + suffix

    def _init_bipartite(self) -> None:
        """Build the bipartite graph (dispatches by mode)."""
        if self.mode == "random":
            self._init_bipartite_random()
        else:
            self._init_bipartite_from_degrees()

    def _init_bipartite_random(self) -> None:
        """Build ER-style random bipartite graph."""
        self.G = bipartite.random_graph(self.n1, self.n2, self.p)
        self._extract_partitions()

    def _init_bipartite_from_degrees(self) -> None:
        """Build bipartite configuration model from degree sequences."""
        self.G = bipartite.configuration_model(
            self.top_degrees,
            self.bottom_degrees,
            create_using=nx.Graph(),
        )
        # Remove parallel edges and self-loops
        self.G = nx.Graph(self.G)
        self._extract_partitions()

    def _extract_partitions(self) -> None:
        """Extract top/bottom node sets from bipartite attribute."""
        self.top_nodes = {
            n
            for n, d in self.G.nodes(data=True)
            if d["bipartite"] == 0
        }
        self.bottom_nodes = {
            n
            for n, d in self.G.nodes(data=True)
            if d["bipartite"] == 1
        }

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        if self.mode == "degree_sequence":
            return f"n1={self.n1}_n2={self.n2}_deg"
        return f"n1={self.n1}_n2={self.n2}_p={self.p:.3g}"

    def get_expected_num_nodes(self) -> int:
        """
        Return the expected number of nodes.

        Returns
        -------
        int
            Total nodes (n1 + n2).
        """
        return self.n1 + self.n2

    def get_expected_num_edges(self) -> int:
        """
        Return the expected number of edges.

        Returns
        -------
        int
            For random mode: ``n1 * n2 * p``.
            For degree-sequence mode: ``sum(top_degrees)``.
        """
        if self.mode == "degree_sequence":
            return sum(self.top_degrees)
        return int(self.n1 * self.n2 * self.p)

    def get_biadjacency_matrix(
        self, sparse: bool = False
    ) -> np.ndarray:
        """
        Return the biadjacency matrix of the bipartite graph.

        The biadjacency matrix B has shape (n1, n2) where B[i,j] = 1
        if there is an edge between top node i and bottom node j.

        Parameters
        ----------
        sparse : bool, default False
            If True, return a scipy sparse matrix.

        Returns
        -------
        ndarray or sparse matrix
            Biadjacency matrix of shape (n1, n2).
        """
        row_order = sorted(self.top_nodes)
        col_order = sorted(self.bottom_nodes)

        if sparse:
            return bipartite.biadjacency_matrix(
                self.G,
                row_order=row_order,
                column_order=col_order,
            )
        else:
            return bipartite.biadjacency_matrix(
                self.G,
                row_order=row_order,
                column_order=col_order,
            ).toarray()

    def get_projection(self, which: str = "top") -> nx.Graph:
        """
        Return a one-mode projection of the bipartite graph.

        The projection connects two nodes if they share a common
        neighbor in the other partition.

        Parameters
        ----------
        which : str, default "top"
            Which partition to project onto: ``"top"`` or ``"bottom"``.

        Returns
        -------
        nx.Graph
            Projected one-mode graph.
        """
        if which == "top":
            return bipartite.projected_graph(
                self.G, self.top_nodes
            )
        elif which == "bottom":
            return bipartite.projected_graph(
                self.G, self.bottom_nodes
            )
        else:
            raise ValueError(
                f"which must be 'top' or 'bottom', got {which}"
            )

    def get_density(self) -> float:
        """
        Return the edge density of the bipartite graph.

        Returns
        -------
        float
            Density = edges / (n1 * n2).
        """
        return self.G.number_of_edges() / (self.n1 * self.n2)

    def is_bipartite(self) -> bool:
        """
        Verify that the graph is bipartite.

        Returns
        -------
        bool
            True if graph is bipartite.
        """
        return bipartite.is_bipartite(self.G)
