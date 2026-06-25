"""
BipartiteGraphGT: Native graph-tool implementation of signed bipartite graphs.

Unified class supporting both ER-style random bipartite graphs and
bipartite configuration model graphs from prescribed degree sequences.

Example
-------
>>> from lrgsglib.graphs.gt import BipartiteGraphGT
>>> # Mode 1: ER random bipartite
>>> bg = BipartiteGraphGT(n1=20, n2=30, p=0.1, pflip=0.1, seed=42)
>>> bg.flip_random_fract_edges()
>>> print(f"Nodes: {bg.N}, Edges: {bg.num_edges}")
>>>
>>> # Mode 2: Configuration model from degree sequences
>>> bg2 = BipartiteGraphGT(
...     top_degrees=[3, 3, 2, 2],
...     bottom_degrees=[2, 2, 2, 2, 2],
...     seed=42,
... )
>>> print(f"Nodes: {bg2.N}, Mode: {bg2.mode}")
"""
from __future__ import annotations

from typing import Optional, Sequence, Set, Tuple

import numpy as np
from numpy.typing import NDArray

try:
    import graph_tool as gt
    from graph_tool import Graph

    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object

from ..SignedGraphGT import SignedGraphGT


__all__ = ["BipartiteGraphGT"]


# Constants
BP_PHTABB = "bp"
BP_SGPATH = ""
BP_STDFN = ""


class BipartiteGraphGT(SignedGraphGT):
    """
    Native graph-tool implementation of signed bipartite (two-mode) graphs.

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
    pflip : float, default 0.0
        Fraction of edges to flip to negative.
    seed : int, optional
        Random seed for reproducibility.
    **kwargs
        Additional arguments (for API compatibility).

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
        Nodes in the first partition (0 to n1-1).
    bottom_nodes : set
        Nodes in the second partition (n1 to n1+n2-1).
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

    The bipartite structure is stored as a vertex property:
    - ``G.vp.bipartite[v]`` is 0 for top nodes, 1 for bottom nodes

    Examples
    --------
    >>> from lrgsglib.graphs.gt import BipartiteGraphGT
    >>> # Random bipartite graph
    >>> bg = BipartiteGraphGT(n1=50, n2=100, p=0.05, seed=42)
    >>> print(f"Top: {len(bg.top_nodes)}, Bottom: {len(bg.bottom_nodes)}")
    Top: 50, Bottom: 100

    >>> # Degree-sequence bipartite graph
    >>> bg2 = BipartiteGraphGT(
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
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ) -> None:
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # --- Resolve build mode ---
        has_deg = top_degrees is not None and bottom_degrees is not None
        has_partial_deg = (
            (top_degrees is not None) != (bottom_degrees is not None)
        )

        if has_partial_deg:
            raise ValueError(
                "Both top_degrees and bottom_degrees must be provided "
                "together for degree-sequence mode."
            )

        if has_deg:
            # Degree-sequence mode
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
            if any(d < 0 for d in self.top_degrees):
                raise ValueError("All top degrees must be non-negative")
            if any(d < 0 for d in self.bottom_degrees):
                raise ValueError(
                    "All bottom degrees must be non-negative"
                )

            self.n1 = len(self.top_degrees)
            self.n2 = len(self.bottom_degrees)
            self.p = 0.0
        else:
            # Random mode
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

        # Store common parameters
        self.only_const_mode = only_const_mode
        self.syshape = (self.n1, self.n2)
        if self.mode == "random":
            self.syshapePth = (
                f"n1={self.n1}_n2={self.n2}_p={self.p:.3g}"
            )
        else:
            self.syshapePth = f"n1={self.n1}_n2={self.n2}_deg"

        # Initialize standard filename
        self._init_std_fname(stdFnameSFFX)
        self.sgpathn = BP_PHTABB if not sgpathn else sgpathn

        # Initialize RNG before generation
        self._rng = np.random.default_rng(seed)

        # Seed graph-tool RNG if seed provided
        if seed is not None:
            gt.seed_rng(seed)

        # Generate the graph
        if not only_const_mode:
            G = self._generate()
        else:
            G = Graph(directed=False)
            sign_prop = G.new_edge_property("int", val=1)
            G.edge_properties["sign"] = sign_prop
            self.top_nodes: Set[int] = set()
            self.bottom_nodes: Set[int] = set()

        # Initialize base class with generated graph
        super().__init__(G=G, pflip=pflip, seed=seed, **kwargs)


    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        if self.mode == "degree_sequence":
            self.std_fname = BP_PHTABB + "_deg" + suffix
        else:
            self.std_fname = BP_PHTABB + suffix

    def _generate(self) -> "Graph":
        """Dispatch to the appropriate generation method."""
        if self.mode == "random":
            return self._generate_random()
        else:
            return self._generate_from_degrees()

    def _generate_random(self) -> "Graph":
        """Generate bipartite graph using ER-style random edges.

        Creates vertices, assigns partition labels, and adds random
        edges between partitions with probability p.

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign and bipartite
            properties.
        """
        G = gt.Graph(directed=False)
        G.add_vertex(self.n1 + self.n2)

        # Add bipartite partition property using numpy array
        bipartite_prop = G.new_vertex_property("int")
        bp_array = np.concatenate([
            np.zeros(self.n1, dtype=int),
            np.ones(self.n2, dtype=int),
        ])
        bipartite_prop.a = bp_array
        G.vertex_properties["bipartite"] = bipartite_prop

        # Store node sets
        self.top_nodes = set(range(self.n1))
        self.bottom_nodes = set(range(self.n1, self.n1 + self.n2))

        # Generate random edges using vectorized numpy operations
        n_potential_edges = self.n1 * self.n2
        random_vals = self._rng.random(n_potential_edges)
        mask = random_vals < self.p

        if np.any(mask):
            edge_indices = np.where(mask)[0]
            top_indices = edge_indices // self.n2
            bottom_indices = edge_indices % self.n2 + self.n1
            edge_list = np.column_stack(
                [top_indices, bottom_indices]
            )
            G.add_edge_list(edge_list)

        # Add sign property (all +1 initially)
        sign_prop = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign_prop

        return G

    def _generate_from_degrees(self) -> "Graph":
        """Generate bipartite configuration model.

        Uses stub-matching approach: create stubs for each node
        according to its prescribed degree, then randomly match top
        stubs with bottom stubs.

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign and bipartite
            properties.
        """
        G = gt.Graph(directed=False)
        G.add_vertex(self.n1 + self.n2)

        # Add bipartite partition property
        bipartite_prop = G.new_vertex_property("int")
        for v in range(self.n1):
            bipartite_prop[v] = 0
        for v in range(self.n1, self.n1 + self.n2):
            bipartite_prop[v] = 1
        G.vertex_properties["bipartite"] = bipartite_prop

        # Store node sets
        self.top_nodes = set(range(self.n1))
        self.bottom_nodes = set(range(self.n1, self.n1 + self.n2))

        # Create stubs for configuration model
        top_stubs = []
        for v, deg in enumerate(self.top_degrees):
            top_stubs.extend([v] * deg)

        bottom_stubs = []
        for v, deg in enumerate(self.bottom_degrees):
            bottom_stubs.extend([v + self.n1] * deg)

        # Shuffle stubs
        self._rng.shuffle(top_stubs)
        self._rng.shuffle(bottom_stubs)

        # Match stubs to create edges (set avoids multi-edges)
        edge_set = set()
        for t_stub, b_stub in zip(top_stubs, bottom_stubs):
            edge_set.add((t_stub, b_stub))

        edge_list = list(edge_set)
        if edge_list:
            G.add_edge_list(edge_list)

        # Add sign property (all +1 initially)
        sign_prop = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign_prop

        return G

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

    @property
    def N(self) -> int:
        """Number of nodes."""
        return self.G.num_vertices()

    def get_actual_top_degrees(self) -> NDArray[np.int_]:
        """
        Return the actual degree sequence of top nodes.

        Only meaningful in degree-sequence mode.

        Returns
        -------
        np.ndarray
            Array of actual degrees for top nodes.

        Raises
        ------
        AttributeError
            If called in random mode.
        """
        if self.mode != "degree_sequence":
            raise AttributeError(
                "get_actual_top_degrees() is only available in "
                "degree-sequence mode."
            )
        return np.array([
            self.G.vertex(v).out_degree() for v in range(self.n1)
        ])

    def get_actual_bottom_degrees(self) -> NDArray[np.int_]:
        """
        Return the actual degree sequence of bottom nodes.

        Only meaningful in degree-sequence mode.

        Returns
        -------
        np.ndarray
            Array of actual degrees for bottom nodes.

        Raises
        ------
        AttributeError
            If called in random mode.
        """
        if self.mode != "degree_sequence":
            raise AttributeError(
                "get_actual_bottom_degrees() is only available in "
                "degree-sequence mode."
            )
        return np.array([
            self.G.vertex(v).out_degree()
            for v in range(self.n1, self.n1 + self.n2)
        ])

    def verify_degrees(self, tolerance: float = 0.1) -> bool:
        """
        Verify that actual degrees match prescribed degrees.

        Due to multi-edge removal, actual degrees may be less than
        prescribed. Only meaningful in degree-sequence mode.

        Parameters
        ----------
        tolerance : float, default 0.1
            Maximum allowed fraction difference.

        Returns
        -------
        bool
            True if degrees match within tolerance.

        Raises
        ------
        AttributeError
            If called in random mode.
        """
        if self.mode != "degree_sequence":
            raise AttributeError(
                "verify_degrees() is only available in "
                "degree-sequence mode."
            )
        actual_top = self.get_actual_top_degrees()
        actual_bottom = self.get_actual_bottom_degrees()
        prescribed_top = np.array(self.top_degrees)
        prescribed_bottom = np.array(self.bottom_degrees)

        top_diff = (
            np.abs(actual_top - prescribed_top).sum()
            / prescribed_top.sum()
        )
        bottom_diff = (
            np.abs(actual_bottom - prescribed_bottom).sum()
            / prescribed_bottom.sum()
        )

        return top_diff <= tolerance and bottom_diff <= tolerance

    def get_biadjacency_matrix(
        self, sparse: bool = False
    ) -> NDArray[np.int_]:
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
        B = np.zeros((self.n1, self.n2), dtype=int)

        for e in self.G.edges():
            i, j = int(e.source()), int(e.target())
            if i >= self.n1:
                i, j = j, i
            B[i, j - self.n1] = 1

        if sparse:
            from scipy import sparse as sp

            return sp.csr_matrix(B)

        return B

    def get_projection(self, which: str = "top") -> "Graph":
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
        Graph
            Projected one-mode graph-tool Graph.
        """
        B = self.get_biadjacency_matrix()

        if which == "top":
            adj = B @ B.T
            n = self.n1
        elif which == "bottom":
            adj = B.T @ B
            n = self.n2
        else:
            raise ValueError(
                f"which must be 'top' or 'bottom', got {which}"
            )

        G_proj = gt.Graph(directed=False)
        G_proj.add_vertex(n)

        edge_list = []
        for i in range(n):
            for j in range(i + 1, n):
                if adj[i, j] > 0:
                    edge_list.append((i, j))

        if edge_list:
            G_proj.add_edge_list(edge_list)

        return G_proj

    def get_density(self) -> float:
        """
        Return the edge density of the bipartite graph.

        Returns
        -------
        float
            Density = edges / (n1 * n2).
        """
        return self.G.num_edges() / (self.n1 * self.n2)

    def is_bipartite(self) -> bool:
        """
        Verify that the graph is bipartite.

        Returns
        -------
        bool
            True if graph is bipartite (no within-partition edges).
        """
        for e in self.G.edges():
            i, j = int(e.source()), int(e.target())
            bp_i = self.G.vp.bipartite[i]
            bp_j = self.G.vp.bipartite[j]
            if bp_i == bp_j:
                return False
        return True

    def __repr__(self) -> str:
        if self.mode == "degree_sequence":
            return (
                f"BipartiteGraphGT(n1={self.n1}, n2={self.n2}, "
                f"mode='degree_sequence', "
                f"edges={self.G.num_edges()})"
            )
        return (
            f"BipartiteGraphGT(n1={self.n1}, n2={self.n2}, "
            f"p={self.p:.3g}, edges={self.G.num_edges()})"
        )

    def __str__(self) -> str:
        return self.__repr__()
