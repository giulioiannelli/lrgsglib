"""
MultiplicativeCascadeGraphGT - Native graph-tool implementation of multiplicative cascade graphs.

This module provides a native graph-tool implementation that constructs multiplicative
cascade graphs directly without falling back to NetworkX.
"""

from __future__ import annotations

from typing import Optional

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

from ..MultispectralGraphGT import MultispectralGraphGT
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


class MultiplicativeCascadeGraphGT(MultispectralGraphGT):
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
    >>> from lrgsglib.graphs.gt import MultiplicativeCascadeGraphGT
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
        super().__init__(G=G, pflip=pflip, seed=seed, sgpathn=sgpathn, **kwargs)


    def _generate(self) -> "Graph":
        """Generate multiplicative cascade graph using native graph-tool operations.

        Steps:
        1) Compute probability matrix via Kronecker products
        2) Select nodes using exponential clocks or standard algorithm
        3) Build edges via vectorized numpy grid operations
        4) Attach grid positions as vertex property
        5) Extract giant component
        6) Add sign property

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign and pos properties.
        """
        # Step 1: Compute probability matrix (uses numpy, same as NX)
        old_state = np.random.get_state()
        np.random.seed(self._rng.integers(0, 2**31))
        self.probability_matrix = multiplicative_cascade_probability_matrix(
            self.p1, self.p2, self.p3, self.p4,
            iterations=self.iterations,
            stochastic=self.stochastic,
        )
        np.random.set_state(old_state)

        size = int(self.probability_matrix.shape[0])

        # Step 2: Select nodes
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

        # Step 3: Build edges (C++ backend if available, else vectorized numpy)
        n_selected = len(selected_indices)
        G = gt.Graph(directed=False)
        G.add_vertex(n_selected)

        try:
            from .cpp import build_cascade_edges
            build_cascade_edges(
                G,
                selected_indices[:, 0].tolist(),
                selected_indices[:, 1].tolist(),
                n_selected,
                size,
                self.periodic,
            )
        except ImportError:
            edge_array = self._build_edges_vectorized(selected_indices, size)
            if len(edge_array) > 0:
                G.add_edge_list(edge_array)

        # Step 4: Attach grid positions as vertex property
        pos = G.new_vertex_property("vector<double>")
        rows = selected_indices[:, 0].astype(np.float64)
        cols = selected_indices[:, 1].astype(np.float64)
        for v_idx in range(n_selected):
            pos[v_idx] = [cols[v_idx], -rows[v_idx]]
        G.vertex_properties["pos"] = pos

        # Step 5: Extract giant component
        if G.num_edges() > 0:
            G = self._extract_giant_component(G)

        if G.num_vertices() == 0:
            raise ValueError(
                "The multiplicative cascade graph is empty after "
                "giant component extraction."
            )

        # Step 6: Add sign property
        sign_prop = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign_prop

        return G

    def _build_edges_vectorized(
        self,
        selected_indices: NDArray[np.integer],
        size: int,
    ) -> NDArray[np.integer]:
        """Build grid edges between selected nodes using vectorized numpy.

        Uses boolean occupancy grids and ``np.roll`` / slicing to find
        right-neighbor and down-neighbor edges in O(size^2) time instead
        of O(n_selected) dict-lookup time.

        Parameters
        ----------
        selected_indices : np.ndarray
            Shape ``(n_selected, 2)`` with ``(row, col)`` of each node.
        size : int
            Grid side length.

        Returns
        -------
        np.ndarray
            Shape ``(n_edges, 2)`` edge list in vertex-index space.
        """
        n_selected = len(selected_indices)

        # Build occupancy grid and index grid
        # grid[i, j] = True if (i, j) is a selected node
        grid = np.zeros((size, size), dtype=bool)
        idx_grid = np.full((size, size), -1, dtype=np.int64)

        rows = selected_indices[:, 0]
        cols = selected_indices[:, 1]
        grid[rows, cols] = True
        idx_grid[rows, cols] = np.arange(n_selected)

        edges = []

        # Right-neighbor edges (horizontal: j -> j+1)
        if self.periodic:
            right_grid = np.roll(grid, -1, axis=1)
            right_idx = np.roll(idx_grid, -1, axis=1)
        else:
            # Non-periodic: only columns 0..size-2 can have a right neighbor
            right_grid = np.zeros_like(grid)
            right_idx = np.full_like(idx_grid, -1)
            right_grid[:, :-1] = grid[:, 1:]
            right_idx[:, :-1] = idx_grid[:, 1:]

        h_mask = grid & right_grid
        if h_mask.any():
            h_rows, h_cols = np.nonzero(h_mask)
            src = idx_grid[h_rows, h_cols]
            dst = right_idx[h_rows, h_cols]
            edges.append(np.column_stack((src, dst)))

        # Down-neighbor edges (vertical: i -> i+1)
        if self.periodic:
            down_grid = np.roll(grid, -1, axis=0)
            down_idx = np.roll(idx_grid, -1, axis=0)
        else:
            down_grid = np.zeros_like(grid)
            down_idx = np.full_like(idx_grid, -1)
            down_grid[:-1, :] = grid[1:, :]
            down_idx[:-1, :] = idx_grid[1:, :]

        v_mask = grid & down_grid
        if v_mask.any():
            v_rows, v_cols = np.nonzero(v_mask)
            src = idx_grid[v_rows, v_cols]
            dst = down_idx[v_rows, v_cols]
            edges.append(np.column_stack((src, dst)))

        if edges:
            return np.vstack(edges)
        return np.empty((0, 2), dtype=np.int64)

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

        Uses vectorized PropertyMap ``.a`` accessor instead of
        per-vertex Python iteration.

        Parameters
        ----------
        G : Graph
            Input graph (may have 'pos' vertex property attached).

        Returns
        -------
        Graph
            Graph containing only the giant component. Internal
            vertex properties (including 'pos') survive pruning.
        """
        comp_labels, hist = label_components(G)

        if len(hist) == 0:
            return G

        # Vectorized: set filter via numpy array
        giant_label = int(np.argmax(hist))
        vfilt = G.new_vertex_property("bool")
        vfilt.a = comp_labels.a == giant_label

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
