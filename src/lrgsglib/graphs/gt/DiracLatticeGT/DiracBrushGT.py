"""
DiracBrushGraphGT - Native graph-tool implementation of Dirac brush graphs.

This module provides a native graph-tool implementation that constructs Dirac brush
graphs directly without falling back to NetworkX.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

try:
    import graph_tool as gt
    from graph_tool import Graph
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object

from ....config.const import DBRUSH_STDFN, DBRUSH_SGPATH
from .DiracLatticeGT import DiracLatticeGraphGT


__all__ = ["DiracBrushGraphGT", "DiracBrushGraph"]


class DiracBrushGraphGT(DiracLatticeGraphGT):
    """Dirac brush: 2D base network with 1D fiber networks - graph-tool backend.

    This is a native graph-tool implementation that constructs the graph directly
    using GT's vertex/edge operations, providing better performance for large graphs.

    A Dirac brush is a special case of the Dirac comb where the base is a 2D
    plane/torus and fibers are 1D lines/rings attached to each base node.

    Parameters
    ----------
    base_x : int
        X dimension of the 2D base graph
    base_y : int
        Y dimension of the 2D base graph
    fiber_nodes : int
        Number of nodes in each fiber graph
    periodic : bool
        Use periodic boundary conditions
    pflip : float
        Fraction of edges to flip to negative
    seed : int, optional
        Random seed for reproducibility

    Attributes
    ----------
    base_x : int
        X dimension of base
    base_y : int
        Y dimension of base
    fiber_nodes : int
        Number of nodes per fiber

    Example
    -------
    >>> from lrgsglib.graphs.gt.dirac import DiracBrushGraphGT
    >>> db = DiracBrushGraphGT(base_x=5, base_y=5, fiber_nodes=3, pflip=0.1, seed=42)
    >>> print(f"Nodes: {db.N}, Edges: {db.num_edges}")
    """

    def __init__(
        self,
        base_x: int,
        base_y: int,
        fiber_nodes: int,
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        stdFnameSFFX: str = DBRUSH_STDFN,
        sgpathn: str = DBRUSH_SGPATH,
        **kwargs,
    ):
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # Store parameters
        self.base_x = base_x
        self.base_y = base_y
        self.fiber_nodes = fiber_nodes

        # Path attributes (for compatibility)
        self.stdFnameSFFX = stdFnameSFFX
        self.sgpathn = sgpathn

        # Initialize base class (which calls _generate)
        super().__init__(
            periodic=periodic,
            pflip=pflip,
            seed=seed,
            **kwargs,
        )

    def _generate(self) -> "Graph":
        """Generate Dirac brush graph using native graph-tool operations.

        Structure:
        - Vertices 0 to (base_x * base_y - 1): base nodes
        - Vertices (base_x * base_y) to end: fiber nodes

        Vertex layout:
        - Base node at (i, j) is at vertex index: i * base_y + j
        - Fiber nodes for base (i, j) are at indices:
          base_total + (i * base_y + j) * fiber_nodes to
          base_total + (i * base_y + j + 1) * fiber_nodes - 1

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign property.
        """
        # Validate parameters
        if self.base_x <= 0:
            raise ValueError("base_x must be a positive integer.")
        if self.base_y <= 0:
            raise ValueError("base_y must be a positive integer.")
        if self.fiber_nodes <= 0:
            raise ValueError("fiber_nodes must be a positive integer.")

        # Calculate total nodes
        base_total = self.base_x * self.base_y
        total_nodes = base_total * (1 + self.fiber_nodes)

        # Create graph-tool graph
        G = gt.Graph(directed=False)
        G.add_vertex(total_nodes)

        # Collect all edges
        edge_list = []

        # Step 1: Add base edges (2D grid/torus)
        base_edges = self._create_2d_edges(
            nx=self.base_x,
            ny=self.base_y,
            start_vertex=0,
            periodic=self.periodic,
        )
        edge_list.extend(base_edges)

        # Step 2: Add fiber edges and anchor connections
        for i in range(self.base_x):
            for j in range(self.base_y):
                # Base vertex index
                base_idx = i * self.base_y + j

                # Fiber starting vertex index
                fiber_start = base_total + base_idx * self.fiber_nodes

                # Add fiber edges (1D line/ring)
                fiber_edges = self._create_1d_edges(
                    n_nodes=self.fiber_nodes,
                    start_vertex=fiber_start,
                    periodic=self.periodic,
                )
                edge_list.extend(fiber_edges)

                # Add anchor edge: connect fiber[0] to base node
                anchor_vertex = fiber_start  # fiber_idx=0
                edge_list.append((base_idx, anchor_vertex))

        # Batch add edges (more efficient)
        if edge_list:
            G.add_edge_list(edge_list)

        # Add sign property (all +1 initially)
        sign_prop = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign_prop

        # Store metadata
        self.dirac_structure = {
            'base_x': self.base_x,
            'base_y': self.base_y,
            'base_nodes': base_total,
            'fiber_nodes': self.fiber_nodes,
            'total_nodes': total_nodes,
            'structure': 'dirac_brush',
            'periodic': self.periodic,
        }

        return G

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        # Total = base nodes + (base_nodes * fiber_nodes)
        base_nodes = self.base_x * self.base_y
        return base_nodes * (1 + self.fiber_nodes)

    def get_base_vertex_index(self, i: int, j: int) -> int:
        """Return vertex index for base node at (i, j).

        Parameters
        ----------
        i : int
            X coordinate in base
        j : int
            Y coordinate in base

        Returns
        -------
        int
            Vertex index
        """
        if i < 0 or i >= self.base_x:
            raise ValueError(f"i must be in [0, {self.base_x})")
        if j < 0 or j >= self.base_y:
            raise ValueError(f"j must be in [0, {self.base_y})")
        return i * self.base_y + j

    def get_base_vertex_indices(self) -> list:
        """Return list of all base vertex indices."""
        base_total = self.base_x * self.base_y
        return list(range(base_total))

    def get_fiber_vertex_indices(self, i: int, j: int) -> list:
        """Return list of fiber vertex indices for base node at (i, j).

        Parameters
        ----------
        i : int
            X coordinate of the base node
        j : int
            Y coordinate of the base node

        Returns
        -------
        list
            List of vertex indices in the fiber attached to this base node
        """
        if i < 0 or i >= self.base_x:
            raise ValueError(f"i must be in [0, {self.base_x})")
        if j < 0 or j >= self.base_y:
            raise ValueError(f"j must be in [0, {self.base_y})")

        base_total = self.base_x * self.base_y
        base_idx = i * self.base_y + j
        fiber_start = base_total + base_idx * self.fiber_nodes
        return list(range(fiber_start, fiber_start + self.fiber_nodes))

    def __repr__(self) -> str:
        return (
            f"DiracBrushGraphGT(base={self.base_x}x{self.base_y}, "
            f"fiber={self.fiber_nodes}, N={self.N})"
        )


# Backward compatibility alias
DiracBrushGraph = DiracBrushGraphGT
