"""
SignedGraphGT: graph-tool backend for signed graphs.

This module provides a graph-tool implementation of SignedGraph that
mirrors the nx_patches.SignedGraph API while leveraging graph-tool's
C++ performance.

Example
-------
>>> from lrgsglib.gt_patches import SignedGraphGT
>>> sg = SignedGraphGT(n=100)  # Create from scratch
>>> # Or convert from NX
>>> from lrgsglib.gt_patches.converters import nx_to_gt
>>> sg_gt = SignedGraphGT.from_networkx(nx_signed_graph)
"""

from typing import Any, Optional, Union

import numpy as np

try:
    import graph_tool as gt
    from graph_tool import Graph
    from graph_tool.spectral import adjacency, laplacian
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object  # Placeholder for type hints

import networkx as nx

from .converters import nx_to_gt, gt_to_nx, get_laplacian_matrix_gt


class SignedGraphGT:
    """
    Signed graph using graph-tool backend.

    This class provides a graph-tool based implementation of signed graphs,
    offering significant performance improvements for large graphs while
    maintaining API compatibility with nx_patches.SignedGraph.

    Parameters
    ----------
    G : Graph, optional
        Existing graph-tool Graph to wrap.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping.
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    G : Graph
        The underlying graph-tool Graph.
    N : int
        Number of nodes.
    pflip : float
        Fraction of edges marked for flipping.

    Notes
    -----
    Key differences from nx_patches.SignedGraph:
    - Uses graph-tool property maps instead of NetworkX edge attributes
    - Spectral operations are faster via graph-tool's C++ backend
    - Memory-efficient for large graphs

    Examples
    --------
    >>> from lrgsglib.gt_patches import SignedGraphGT
    >>> # Create and flip edges
    >>> sg = SignedGraphGT.from_complete(n=50, pflip=0.2, seed=42)
    >>> sg.flip_random_fract_edges()
    >>> print(f"Negative edges: {sg.count_negative_edges()}")
    """

    def __init__(
        self,
        G: Optional["Graph"] = None,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> None:
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        self.pflip = pflip
        self.seed = seed
        self._rng = np.random.default_rng(seed)

        if G is not None:
            self.G = G
        else:
            self.G = Graph(directed=False)

        # Ensure sign property exists
        if "sign" not in self.G.edge_properties:
            self.G.edge_properties["sign"] = self.G.new_edge_property("int", val=1)

        # Track edges marked for flipping
        self._flip_edges = set()

    @property
    def N(self) -> int:
        """Number of nodes."""
        return self.G.num_vertices()

    @property
    def num_edges(self) -> int:
        """Number of edges."""
        return self.G.num_edges()

    @classmethod
    def from_networkx(
        cls,
        G_nx: nx.Graph,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> "SignedGraphGT":
        """
        Create SignedGraphGT from NetworkX graph.

        Parameters
        ----------
        G_nx : nx.Graph
            NetworkX graph (with optional 'sign' edge attribute).
        pflip : float, default 0.0
            Fraction of edges to mark for flipping.
        seed : int, optional
            Random seed.

        Returns
        -------
        SignedGraphGT
            New instance wrapping converted graph.
        """
        G_gt = nx_to_gt(G_nx)
        return cls(G=G_gt, pflip=pflip, seed=seed)

    @classmethod
    def from_complete(
        cls,
        n: int,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> "SignedGraphGT":
        """
        Create complete (fully connected) signed graph.

        Parameters
        ----------
        n : int
            Number of nodes.
        pflip : float, default 0.0
            Fraction of edges to mark for flipping.
        seed : int, optional
            Random seed.

        Returns
        -------
        SignedGraphGT
            Complete graph with all positive edges.
        """
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")

        from graph_tool.generation import complete_graph
        G = complete_graph(n)
        return cls(G=G, pflip=pflip, seed=seed)

    @classmethod
    def from_lattice(
        cls,
        shape: tuple,
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> "SignedGraphGT":
        """
        Create lattice signed graph.

        Parameters
        ----------
        shape : tuple
            Shape of lattice (e.g., (10, 10) for 2D).
        periodic : bool, default True
            Whether to use periodic boundary conditions.
        pflip : float, default 0.0
            Fraction of edges to mark for flipping.
        seed : int, optional
            Random seed.

        Returns
        -------
        SignedGraphGT
            Lattice graph with all positive edges.
        """
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")

        from graph_tool.generation import lattice
        G = lattice(shape, periodic=periodic)
        return cls(G=G, pflip=pflip, seed=seed)

    def to_networkx(self) -> nx.Graph:
        """
        Convert to NetworkX graph.

        Returns
        -------
        nx.Graph
            NetworkX graph with 'sign' edge attribute.
        """
        return gt_to_nx(self.G)

    def mark_flip_edges(self, fraction: Optional[float] = None) -> None:
        """
        Mark a fraction of edges for sign flipping.

        Parameters
        ----------
        fraction : float, optional
            Fraction to mark. If None, uses self.pflip.
        """
        frac = fraction if fraction is not None else self.pflip
        if frac <= 0:
            return

        edges = list(self.G.edges())
        n_flip = int(len(edges) * frac)

        if n_flip > 0:
            indices = self._rng.choice(len(edges), size=n_flip, replace=False)
            self._flip_edges = {edges[i] for i in indices}

    def flip_random_fract_edges(self, fraction: Optional[float] = None) -> None:
        """
        Flip signs of a random fraction of edges.

        Parameters
        ----------
        fraction : float, optional
            Fraction to flip. If None, uses self.pflip.
        """
        self.mark_flip_edges(fraction)
        sign_prop = self.G.edge_properties["sign"]

        for e in self._flip_edges:
            sign_prop[e] = -1

    def flip_edge(self, u: int, v: int) -> None:
        """
        Flip sign of a specific edge.

        Parameters
        ----------
        u, v : int
            Vertex indices.
        """
        e = self.G.edge(u, v)
        if e is not None:
            self.G.edge_properties["sign"][e] *= -1

    def get_edge_sign(self, u: int, v: int) -> int:
        """
        Get sign of an edge.

        Parameters
        ----------
        u, v : int
            Vertex indices.

        Returns
        -------
        int
            Edge sign (+1 or -1), or 0 if edge doesn't exist.
        """
        e = self.G.edge(u, v)
        if e is not None:
            return self.G.edge_properties["sign"][e]
        return 0

    def count_negative_edges(self) -> int:
        """
        Count number of negative edges.

        Returns
        -------
        int
            Number of edges with sign = -1.
        """
        sign_prop = self.G.edge_properties["sign"]
        return sum(1 for e in self.G.edges() if sign_prop[e] == -1)

    def count_positive_edges(self) -> int:
        """
        Count number of positive edges.

        Returns
        -------
        int
            Number of edges with sign = +1.
        """
        return self.num_edges - self.count_negative_edges()

    def get_signed_adjacency(self) -> np.ndarray:
        """
        Get signed adjacency matrix.

        Returns
        -------
        np.ndarray
            Signed adjacency matrix A where A[i,j] = sign of edge (i,j).
        """
        return adjacency(self.G, weight=self.G.edge_properties["sign"]).toarray()

    def get_signed_laplacian(self) -> np.ndarray:
        """
        Get signed Laplacian matrix.

        Returns
        -------
        np.ndarray
            Signed Laplacian L = D - A_signed.
        """
        return get_laplacian_matrix_gt(self.G, signed=True)

    def get_laplacian_spectrum(self, k: Optional[int] = None) -> np.ndarray:
        """
        Compute Laplacian eigenvalues.

        Parameters
        ----------
        k : int, optional
            Number of eigenvalues to compute. If None, compute all.

        Returns
        -------
        np.ndarray
            Sorted eigenvalues.
        """
        L = self.get_signed_laplacian()

        if k is None or k >= self.N:
            # Full spectrum
            eigenvalues = np.linalg.eigvalsh(L)
        else:
            # Partial spectrum using scipy sparse
            from scipy.sparse.linalg import eigsh
            from scipy.sparse import csr_matrix

            L_sparse = csr_matrix(L)
            eigenvalues, _ = eigsh(L_sparse, k=k, which='SM')

        return np.sort(eigenvalues)

    def get_frustration_index(self) -> float:
        """
        Compute frustration index (fraction of frustrated triangles).

        Returns
        -------
        float
            Fraction of triangles that are frustrated (odd # negative edges).
        """
        # Count triangles and frustrated triangles
        # A triangle is frustrated if it has 1 or 3 negative edges
        sign_prop = self.G.edge_properties["sign"]
        n_triangles = 0
        n_frustrated = 0

        for v in self.G.vertices():
            neighbors = list(v.out_neighbors())
            for i, n1 in enumerate(neighbors):
                for n2 in neighbors[i + 1:]:
                    # Check if n1-n2 edge exists (completing triangle)
                    e12 = self.G.edge(n1, n2)
                    if e12 is not None:
                        n_triangles += 1
                        # Count negative edges in triangle
                        e_v_n1 = self.G.edge(v, n1)
                        e_v_n2 = self.G.edge(v, n2)
                        neg_count = sum([
                            sign_prop[e_v_n1] == -1,
                            sign_prop[e_v_n2] == -1,
                            sign_prop[e12] == -1,
                        ])
                        if neg_count % 2 == 1:  # 1 or 3 negative edges
                            n_frustrated += 1

        # Each triangle counted 3 times (once per vertex)
        n_triangles //= 3
        n_frustrated //= 3

        return n_frustrated / n_triangles if n_triangles > 0 else 0.0

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"SignedGraphGT(N={self.N}, edges={self.num_edges}, "
            f"negative={neg}, pflip={self.pflip})"
        )
