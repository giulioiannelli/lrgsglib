"""
BipartiteGraphNX: Signed bipartite (two-mode) network.

Bipartite graphs have two disjoint sets of nodes where edges only connect
nodes from different sets. Common applications include:
- Author-paper networks
- User-item recommendation systems
- Protein-drug interaction networks

Example
-------
>>> from lrgsglib.graphs.nx import BipartiteGraphNX
>>> bg = BipartiteGraphNX(n1=20, n2=30, p=0.1, pflip=0.1, seed=42)
>>> bg.flip_random_fract_edges()
>>> print(f"Nodes: {bg.N}, Edges: {bg.G.number_of_edges()}")
"""

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
    Signed bipartite (two-mode) random graph.

    Creates a bipartite graph with two disjoint node sets where edges
    only connect nodes from different sets. Edge probability can be
    uniform or specified per node pair.

    Parameters
    ----------
    n1 : int
        Number of nodes in the first partition (top nodes).
    n2 : int
        Number of nodes in the second partition (bottom nodes).
    p : float, default 0.1
        Probability of edge between any top-bottom node pair.
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
    n1 : int
        Size of first partition.
    n2 : int
        Size of second partition.
    p : float
        Edge probability.
    top_nodes : set
        Nodes in the first partition.
    bottom_nodes : set
        Nodes in the second partition.

    Notes
    -----
    Bipartite graphs have special spectral properties:
    - Eigenvalues are symmetric around zero
    - The spectrum encodes the singular values of the biadjacency matrix
    - Useful for studying two-mode networks and recommendation systems

    The bipartite structure is stored as node attributes:
    - `G.nodes[v]['bipartite']` is 0 for top nodes, 1 for bottom nodes

    Examples
    --------
    >>> from lrgsglib.graphs.nx import BipartiteGraphNX
    >>> # Create random bipartite graph
    >>> bg = BipartiteGraphNX(n1=50, n2=100, p=0.05, seed=42)
    >>> print(f"Top: {len(bg.top_nodes)}, Bottom: {len(bg.bottom_nodes)}")
    Top: 50, Bottom: 100

    >>> # Access biadjacency matrix
    >>> B = bg.get_biadjacency_matrix()
    >>> print(f"Biadjacency shape: {B.shape}")
    Biadjacency shape: (50, 100)
    """

    def __init__(
        self,
        n1: int,
        n2: int,
        p: float = 0.1,
        sgpathn: str = BP_SGPATH,
        stdFnameSFFX: str = BP_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if n1 < 1 or n2 < 1:
            raise ValueError(f"n1 and n2 must be >= 1, got n1={n1}, n2={n2}")
        if not 0 <= p <= 1:
            raise ValueError(f"p must be in [0, 1], got {p}")

        self.n1 = n1
        self.n2 = n2
        self.p = p
        self.only_const_mode = only_const_mode
        self._init_std_fname(stdFnameSFFX)
        self.sgpathn = BP_PHTABB if not sgpathn else sgpathn
        self.syshape = (n1, n2)
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
        self.std_fname = BP_PHTABB + suffix

    def _init_bipartite(self) -> None:
        """Build the bipartite graph."""
        # Generate random bipartite graph
        self.G = bipartite.random_graph(self.n1, self.n2, self.p)

        # Store node sets
        self.top_nodes = {n for n, d in self.G.nodes(data=True) if d["bipartite"] == 0}
        self.bottom_nodes = {
            n for n, d in self.G.nodes(data=True) if d["bipartite"] == 1
        }

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
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
            Expected edges (n1 * n2 * p).
        """
        return int(self.n1 * self.n2 * self.p)

    def get_biadjacency_matrix(self, sparse: bool = False) -> np.ndarray:
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
                self.G, row_order=row_order, column_order=col_order
            )
        else:
            return bipartite.biadjacency_matrix(
                self.G, row_order=row_order, column_order=col_order
            ).toarray()

    def get_projection(self, which: str = "top") -> nx.Graph:
        """
        Return a one-mode projection of the bipartite graph.

        The projection connects two nodes if they share a common neighbor
        in the other partition.

        Parameters
        ----------
        which : str, default "top"
            Which partition to project onto: "top" or "bottom".

        Returns
        -------
        nx.Graph
            Projected one-mode graph.
        """
        if which == "top":
            return bipartite.projected_graph(self.G, self.top_nodes)
        elif which == "bottom":
            return bipartite.projected_graph(self.G, self.bottom_nodes)
        else:
            raise ValueError(f"which must be 'top' or 'bottom', got {which}")

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
