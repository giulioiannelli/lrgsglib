"""
VicsekGraphGT - Native graph-tool implementation of Vicsek (Palla-Lovasz-Vicsek) graphs.

This module provides a native graph-tool implementation that constructs Vicsek graphs
directly without falling back to NetworkX.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from numpy.typing import NDArray

try:
    import graph_tool as gt
    from graph_tool import Graph
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object

from ..MultispectralGraphGT import MultispectralGraphGT
from ....config.const import VCK_STDFN, VCK_SGPATH
from ...nx.MultispectralGraphNX.generators_msg import (
    initial_measure,
    link_probabilities,
)


__all__ = ["VicsekGraphGT", "VicsekGraph"]


class VicsekGraphGT(MultispectralGraphGT):
    """Vicsek (Palla-Lovasz-Vicsek) graph via Kronecker product - graph-tool backend.

    This is a native graph-tool implementation that constructs the graph directly
    using GT's vertex/edge operations, providing better performance for large graphs.

    Parameters
    ----------
    N : int
        Number of nodes
    k : int
        Number of Kronecker iterations
    pij : np.ndarray, optional
        Initial measure matrix. If not provided, generated from m and symmetric.
    m : int
        Dimension of initial measure matrix (if pij not provided)
    symmetric : bool
        Whether to symmetrize initial measure (if pij not provided)
    pflip : float
        Fraction of edges to flip to negative
    seed : int, optional
        Random seed for reproducibility

    Attributes
    ----------
    N_target : int
        Requested number of nodes
    k : int
        Number of Kronecker iterations
    probability_matrix : np.ndarray
        Final probability matrix P_k after k iterations
    indices : np.ndarray
        Node position indices in probability space

    Example
    -------
    >>> from lrgsglib.graphs.gt.multispectral import VicsekGraphGT
    >>> v = VicsekGraphGT(N=100, k=3, m=3, pflip=0.1, seed=42)
    >>> print(f"Nodes: {v.N}, Edges: {v.num_edges}")
    """

    def __init__(
        self,
        N: int,
        k: int,
        pij: Optional[NDArray[np.floating]] = None,
        m: int = 3,
        symmetric: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        stdFnameSFFX: str = VCK_STDFN,
        sgpathn: str = VCK_SGPATH,
        **kwargs,
    ):
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # Store parameters
        self.N_target = N
        self.k = k
        self.pij_input = pij
        self.m = m
        self.symmetric = symmetric

        # Computed attributes (set during generation)
        self.probability_matrix: Optional[NDArray[np.floating]] = None
        self.indices: Optional[NDArray[np.integer]] = None
        self.sample_fraction = 1.0

        # Path attributes (for compatibility)
        self.stdFnameSFFX = stdFnameSFFX
        self.sgpathn = sgpathn

        # Initialize random state before generation
        self._rng = np.random.default_rng(seed)

        # Generate the graph
        G = self._generate()

        # Initialize base class with generated graph
        super().__init__(G=G, pflip=pflip, seed=seed, sgpathn=sgpathn, **kwargs)


    def _generate(self) -> "Graph":
        """Generate Vicsek graph using native graph-tool operations.

        Steps:
        1) Build the link probability matrix P_k by iterating p_ij k times.
        2) Distribute N nodes uniformly on [0,1].
        3) Connect each pair with probability given by P_k at their coordinates.

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign property.
        """
        # Step 1: Generate or use provided initial measure
        if self.pij_input is None:
            # Use stored RNG for reproducibility
            old_state = np.random.get_state()
            np.random.seed(self._rng.integers(0, 2**31))
            pij = initial_measure(self.m, symmetric=self.symmetric)
            np.random.set_state(old_state)
        else:
            pij = np.asarray(self.pij_input, dtype=float)

        # Step 2: Compute probability matrix via Kronecker products
        Pk = link_probabilities(pij, self.k)
        self.probability_matrix = Pk

        m = pij.shape[0]
        M = int(m ** int(self.k))

        # Step 3: Distribute N nodes uniformly and map to indices
        coords = self._rng.random(int(self.N_target))
        indices = np.floor(coords * M).astype(int)
        indices = np.clip(indices, 0, M - 1)
        self.indices = indices

        # Step 4: Create graph-tool graph
        G = gt.Graph(directed=False)
        G.add_vertex(int(self.N_target))

        # Step 5: Add edges based on probability matrix
        # Pre-generate random values for efficiency
        n_pairs = (self.N_target * (self.N_target - 1)) // 2
        rand_vals = self._rng.random(n_pairs)

        # Add edges where random value < probability
        edge_list = []
        rand_idx = 0
        for i in range(int(self.N_target)):
            for j in range(i + 1, int(self.N_target)):
                p_link = float(Pk[indices[i], indices[j]])
                if rand_vals[rand_idx] < p_link:
                    edge_list.append((i, j))
                rand_idx += 1

        # Batch add edges (more efficient)
        if edge_list:
            G.add_edge_list(edge_list)

        # Step 6: Add sign property (all +1 initially)
        sign_prop = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign_prop

        return G

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        return self.N_target

    @property
    def N(self) -> int:
        """Number of nodes (overrides base to ensure consistency)."""
        return self.G.num_vertices()

    def __repr__(self) -> str:
        return f"VicsekGraphGT(N={self.N}, k={self.k})"


# Backward compatibility alias
VicsekGraph = VicsekGraphGT
