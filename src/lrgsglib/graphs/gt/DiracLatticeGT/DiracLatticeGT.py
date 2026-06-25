"""
DiracLatticeGraphGT - Abstract base class for native graph-tool Dirac lattice structures.

This module provides a native graph-tool implementation of Dirac lattice graphs
that constructs graphs directly without falling back to NetworkX.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

try:
    import graph_tool as gt
    from graph_tool import Graph
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object

from ..SignedGraphGT import SignedGraphGT


__all__ = ["DiracLatticeGraphGT"]


class DiracLatticeGraphGT(SignedGraphGT):
    """Abstract base class for Dirac structure graphs (comb, brush) - graph-tool backend.

    Dirac structures are hierarchical graphs with a base graph and fiber
    graphs attached to each base node. This structure allows efficient
    spectral computation by separating base and fiber spectra.

    This is a native graph-tool implementation that constructs the graph directly
    using GT's vertex/edge operations.

    Parameters
    ----------
    periodic : bool
        Use periodic boundary conditions
    pflip : float
        Fraction of edges to flip to negative
    seed : int, optional
        Random seed for reproducibility

    Attributes
    ----------
    periodic : bool
        Whether periodic boundary conditions are used
    dirac_structure : dict or None
        Metadata about the Dirac structure
    """

    def __init__(
        self,
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ):
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # Store parameters
        self.periodic = periodic

        # Dirac-specific attributes
        self.dirac_structure: Optional[Dict[str, Any]] = None
        self._base_eigenvalues = None
        self._fiber_laplacian = None

        # Initialize random state before generation
        self._rng = np.random.default_rng(seed)

        # Generate the graph (implemented by subclasses)
        G = self._generate()

        # Initialize base class with generated graph
        super().__init__(G=G, pflip=pflip, seed=seed, **kwargs)


    def _generate(self) -> "Graph":
        """Generate the Dirac graph. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement _generate()")

    def _create_1d_edges(
        self,
        n_nodes: int,
        start_vertex: int,
        periodic: bool,
    ) -> list:
        """Create edge list for 1D line/ring structure.

        Parameters
        ----------
        n_nodes : int
            Number of nodes in the line
        start_vertex : int
            Starting vertex index in the global graph
        periodic : bool
            Whether to create a ring (periodic) or line (non-periodic)

        Returns
        -------
        list
            List of (source, target) edge tuples
        """
        edges = []
        for i in range(n_nodes - 1):
            edges.append((start_vertex + i, start_vertex + i + 1))
        if periodic and n_nodes > 2:
            edges.append((start_vertex + n_nodes - 1, start_vertex))
        return edges

    def _create_2d_edges(
        self,
        nx: int,
        ny: int,
        start_vertex: int,
        periodic: bool,
    ) -> list:
        """Create edge list for 2D grid/torus structure.

        Parameters
        ----------
        nx : int
            Number of nodes in x dimension
        ny : int
            Number of nodes in y dimension
        start_vertex : int
            Starting vertex index in the global graph
        periodic : bool
            Whether to create a torus (periodic) or grid (non-periodic)

        Returns
        -------
        list
            List of (source, target) edge tuples
        """
        edges = []

        def idx(i, j):
            return start_vertex + i * ny + j

        # Horizontal edges
        for i in range(nx):
            for j in range(ny - 1):
                edges.append((idx(i, j), idx(i, j + 1)))
            if periodic and ny > 1:
                edges.append((idx(i, ny - 1), idx(i, 0)))

        # Vertical edges
        for i in range(nx - 1):
            for j in range(ny):
                edges.append((idx(i, j), idx(i + 1, j)))
        if periodic and nx > 1:
            for j in range(ny):
                edges.append((idx(nx - 1, j), idx(0, j)))

        return edges

    @property
    def base_eigenvalues(self):
        """Cached base eigenvalues."""
        return self._base_eigenvalues

    @base_eigenvalues.setter
    def base_eigenvalues(self, value):
        self._base_eigenvalues = value

    @property
    def fiber_laplacian(self):
        """Cached fiber Laplacian."""
        return self._fiber_laplacian

    @fiber_laplacian.setter
    def fiber_laplacian(self, value):
        self._fiber_laplacian = value

    def is_dirac_structure(self) -> bool:
        """Always True for DiracLatticeGraphGT subclasses."""
        return True

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement get_expected_num_nodes()")

    @property
    def N(self) -> int:
        """Number of nodes."""
        return self.G.num_vertices()
