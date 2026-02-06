"""
MultiplicativeCascadeGraphGT - Native graph-tool implementation of multiplicative cascade graphs.

This module provides a native graph-tool implementation that constructs multiplicative
cascade graphs directly without falling back to NetworkX.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

try:
    import graph_tool as gt
    from graph_tool import Graph
    from graph_tool.topology import label_components
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object

from ....gt_patches.signed_graph_gt import SignedGraphGT
from ....config.const import (
    MSG_P1,
    MSG_P2,
    MSG_P3,
    MSG_P4,
    MSG_FRACTION,
    MSG_ITERATIONS,
    MC_STDFN,
    MC_SGPATH,
    DEFAULT_P_FSTR_FMT,
)
from ...nx.MultispectralGraphNX.generators_msg import (
    multiplicative_cascade_probability_matrix,
)


__all__ = ["MultiplicativeCascadeGraphGT", "MultiplicativeCascadeGraph"]


class MultiplicativeCascadeGraphGT(SignedGraphGT):
    """Multiplicative cascade graph with hierarchical structure - graph-tool backend.

    This is a native graph-tool implementation that constructs the graph directly
    using GT's vertex/edge operations, providing better performance for large graphs.

    Parameters
    ----------
    p1, p2, p3, p4 : float
        Seed probabilities for cascade quadrants
    fraction : float
        Fraction of nodes to sample
    iterations : int
        Number of cascade iterations
    stochastic : bool
        Use stochastic cascade mode
    periodic : bool
        Use periodic boundary conditions
    variant : str
        Cascade variant: "standard" or "exp_clocks" (default: "exp_clocks")
    pflip : float
        Fraction of edges to flip to negative
    seed : int, optional
        Random seed for reproducibility

    Attributes
    ----------
    probability_matrix : np.ndarray
        The computed probability matrix
    sample_fraction : float
        Fraction of nodes sampled
    seed_probabilities : tuple
        The (p1, p2, p3, p4) seed probabilities

    Example
    -------
    >>> from lrgsglib.graphs.gt.multispectral import MultiplicativeCascadeGraphGT
    >>> mc = MultiplicativeCascadeGraphGT(p1=0.8, p2=0.6, fraction=0.4, iterations=5, seed=42)
    >>> print(f"Nodes: {mc.N}, Edges: {mc.num_edges}")
    """

    def __init__(
        self,
        p1: float = MSG_P1,
        p2: float = MSG_P2,
        p3: float = MSG_P3,
        p4: float = MSG_P4,
        fraction: float = MSG_FRACTION,
        iterations: int = MSG_ITERATIONS,
        stochastic: bool = False,
        periodic: bool = False,
        variant: str = "exp_clocks",
        pflip: float = 0.0,
        seed: Optional[int] = None,
        stdFnameSFFX: str = MC_STDFN,
        sgpathn: str = MC_SGPATH,
        out_suffix: str = "",
        **kwargs,
    ):
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # Store parameters
        self.p1, self.p2, self.p3, self.p4 = p1, p2, p3, p4
        self.fraction = fraction
        self.iterations = iterations
        self.stochastic = stochastic
        self.periodic = periodic
        self.variant = variant
        self.out_suffix = out_suffix

        # Computed attributes
        self.seed_probabilities = (p1, p2, p3, p4)
        self.sample_fraction = fraction
        self.probability_matrix: Optional[NDArray[np.floating]] = None

        # Path attributes (for compatibility)
        self.stdFnameSFFX = stdFnameSFFX
        self.sgpathn = sgpathn

        # Build syshapePth before calling parent
        fmtfunc = lambda x: f"{x:{DEFAULT_P_FSTR_FMT}}"
        mode_str = "stochastic" if stochastic else "regular"
        self.syshapePth = (
            f"i={iterations}_"
            f"f={fmtfunc(fraction)}_"
            f"P={fmtfunc(p1)}_{fmtfunc(p2)}_{fmtfunc(p3)}_{fmtfunc(p4)}_"
            f"{mode_str}"
        )

        # Initialize random state before generation
        self._rng = np.random.default_rng(seed)

        # Generate the graph
        G = self._generate()

        # Initialize base class with generated graph
        super().__init__(G=G, pflip=pflip, seed=seed)

        # Apply sign flips if requested
        if pflip > 0:
            self.flip_random_fract_edges()

    def _generate(self) -> "Graph":
        """Generate multiplicative cascade graph using native graph-tool operations.

        Steps:
        1) Compute probability matrix via Kronecker products
        2) Create 2D grid vertices
        3) Select nodes using exponential clocks algorithm
        4) Build subgraph of selected nodes with grid edges
        5) Extract giant component
        6) Add sign property

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign property.
        """
        # Step 1: Compute probability matrix (uses numpy, same as NX)
        # Use stored RNG for reproducibility
        old_state = np.random.get_state()
        np.random.seed(self._rng.integers(0, 2**31))
        self.probability_matrix = multiplicative_cascade_probability_matrix(
            self.p1, self.p2, self.p3, self.p4,
            iterations=self.iterations,
            stochastic=self.stochastic,
        )
        np.random.set_state(old_state)

        size = int(self.probability_matrix.shape[0])

        # Step 2-3: Select nodes using exp_clocks or standard algorithm
        if self.variant == "exp_clocks":
            selected_indices = self._select_nodes_exp_clocks(
                self.probability_matrix,
                self.fraction,
            )
        else:
            selected_indices = self._select_nodes_standard(
                self.probability_matrix,
                self.fraction,
            )

        if len(selected_indices) == 0:
            raise ValueError("No nodes selected for the cascade graph.")

        # Step 4: Build subgraph with only selected nodes
        # Create mapping from 2D coords to vertex indices
        coord_to_idx = {tuple(coord): i for i, coord in enumerate(selected_indices)}
        n_selected = len(selected_indices)

        # Create graph
        G = gt.Graph(directed=False)
        G.add_vertex(n_selected)

        # Build edge list for grid edges between selected nodes
        edge_list = []
        for idx, (i, j) in enumerate(selected_indices):
            # Check all 4 neighbors (up, down, left, right)
            neighbors = [
                ((i - 1) % size, j) if self.periodic else (i - 1, j),
                ((i + 1) % size, j) if self.periodic else (i + 1, j),
                (i, (j - 1) % size) if self.periodic else (i, j - 1),
                (i, (j + 1) % size) if self.periodic else (i, j + 1),
            ]

            for ni, nj in neighbors:
                # Skip out-of-bounds for non-periodic
                if not self.periodic and (ni < 0 or ni >= size or nj < 0 or nj >= size):
                    continue

                neighbor_coord = (ni, nj)
                if neighbor_coord in coord_to_idx:
                    neighbor_idx = coord_to_idx[neighbor_coord]
                    # Only add edge once (when idx < neighbor_idx)
                    if idx < neighbor_idx:
                        edge_list.append((idx, neighbor_idx))

        if edge_list:
            G.add_edge_list(edge_list)

        # Step 5: Extract giant component
        if G.num_edges() > 0:
            G = self._extract_giant_component(G)

        if G.num_vertices() == 0:
            raise ValueError("The multiplicative cascade graph is empty after giant component extraction.")

        # Step 6: Add sign property
        sign_prop = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign_prop

        return G

    def _select_nodes_exp_clocks(
        self,
        prob_matrix: NDArray[np.floating],
        fraction: float,
    ) -> NDArray[np.integer]:
        """Select nodes using exponential clocks algorithm.

        Parameters
        ----------
        prob_matrix : np.ndarray
            Square probability matrix
        fraction : float
            Fraction of nodes to select

        Returns
        -------
        np.ndarray
            Array of (i, j) coordinates for selected nodes
        """
        N = int(prob_matrix.shape[0])
        k = int(round(float(fraction) * (N ** 2)))
        k = max(1, min(k, N * N))

        # Generate uniform random values
        U = self._rng.random((N, N))

        # Compute exponential clock times
        with np.errstate(divide="ignore", invalid="ignore"):
            time = -np.log1p(-U) / prob_matrix
            time = np.where(prob_matrix > 0.0, time, np.inf)

        # Select k nodes with smallest times
        flat_indices = np.argpartition(time.ravel(), k - 1)[:k]
        indices_2d = np.column_stack(np.unravel_index(flat_indices, (N, N)))

        return indices_2d

    def _select_nodes_standard(
        self,
        prob_matrix: NDArray[np.floating],
        fraction: float,
    ) -> NDArray[np.integer]:
        """Select nodes using standard probability-based sampling.

        Parameters
        ----------
        prob_matrix : np.ndarray
            Square probability matrix
        fraction : float
            Fraction of nodes to select

        Returns
        -------
        np.ndarray
            Array of (i, j) coordinates for selected nodes
        """
        N = int(prob_matrix.shape[0])

        # Sample each cell with probability proportional to prob_matrix * fraction
        thresholds = self._rng.random((N, N))
        selected_mask = thresholds < (prob_matrix * fraction)

        # Get coordinates of selected nodes
        selected_coords = np.argwhere(selected_mask)

        return selected_coords

    def _extract_giant_component(self, G: "Graph") -> "Graph":
        """Extract the giant (largest) connected component.

        Parameters
        ----------
        G : Graph
            Input graph

        Returns
        -------
        Graph
            Graph containing only the giant component
        """
        # Get component labels
        comp_labels, hist = label_components(G)

        if len(hist) == 0:
            return G

        # Find the largest component
        giant_label = int(np.argmax(hist))

        # Create vertex filter for giant component
        vfilt = G.new_vertex_property("bool")
        for v in G.vertices():
            vfilt[v] = (comp_labels[v] == giant_label)

        # Set filter and prune to create new graph
        G.set_vertex_filter(vfilt)
        G_giant = Graph(G, prune=True)

        return G_giant

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        if self.probability_matrix is None:
            return 0
        size = int(self.probability_matrix.shape[0])
        return max(1, int(round(self.sample_fraction * (size ** 2))))

    @property
    def N(self) -> int:
        """Number of nodes."""
        return self.G.num_vertices()

    def __repr__(self) -> str:
        return f"MultiplicativeCascadeGraphGT(N={self.N}, variant={self.variant})"


# Backward compatibility alias
MultiplicativeCascadeGraph = MultiplicativeCascadeGraphGT
