"""
DiracCombGraphGT - Native graph-tool implementation of Dirac comb graphs.

This module provides a native graph-tool implementation that constructs Dirac comb
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

from ....config.const import DCOMB_STDFN, DCOMB_SGPATH
from .DiracLatticeGraphGT import DiracLatticeGraphGT


__all__ = ["DiracCombGraphGT", "DiracCombGraph"]


class DiracCombGraphGT(DiracLatticeGraphGT):
    """Dirac comb: 1D base network with 1D fiber networks - graph-tool backend.

    This is a native graph-tool implementation that constructs the graph directly
    using GT's vertex/edge operations, providing better performance for large graphs.

    A Dirac comb consists of a 1D base graph where each node has a 1D fiber graph
    attached to it. The fiber is connected to the base at one anchor node.

    Parameters
    ----------
    base_nodes : int
        Number of nodes in the base graph
    fiber_nodes : int
        Number of nodes in each fiber graph
    base_type : str
        Type of base graph ('line')
    periodic : bool
        Use periodic boundary conditions
    pflip : float
        Fraction of edges to flip to negative
    seed : int, optional
        Random seed for reproducibility

    Attributes
    ----------
    base_nodes : int
        Number of nodes in base
    fiber_nodes : int
        Number of nodes per fiber
    base_type : str
        Type of base graph

    Example
    -------
    >>> from lrgsglib.graphs.gt.dirac import DiracCombGraphGT
    >>> dc = DiracCombGraphGT(base_nodes=10, fiber_nodes=5, pflip=0.1, seed=42)
    >>> print(f"Nodes: {dc.N}, Edges: {dc.num_edges}")
    """

    def __init__(
        self,
        base_nodes: int,
        fiber_nodes: int,
        base_type: str = "line",
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        stdFnameSFFX: str = DCOMB_STDFN,
        sgpathn: str = DCOMB_SGPATH,
        **kwargs,
    ):
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # Store parameters
        self.base_nodes = base_nodes
        self.fiber_nodes = fiber_nodes
        self.base_type = base_type

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
        """Generate Dirac comb graph using native graph-tool operations.

        Structure:
        - Vertices 0 to base_nodes-1: base nodes
        - Vertices base_nodes to end: fiber nodes (fiber_nodes per base node)

        Vertex layout:
        - Base node i is at vertex index i
        - Fiber nodes for base i are at indices:
          base_nodes + i * fiber_nodes to base_nodes + (i+1) * fiber_nodes - 1

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign property.
        """
        # Validate parameters
        if self.base_nodes <= 0:
            raise ValueError("base_nodes must be a positive integer.")
        if self.fiber_nodes <= 0:
            raise ValueError("fiber_nodes must be a positive integer.")
        if self.base_type != "line":
            raise ValueError(f"Unsupported base_type: {self.base_type}. Use 'line'.")

        # Calculate total nodes
        total_nodes = self.base_nodes * (1 + self.fiber_nodes)

        # Create graph-tool graph
        G = gt.Graph(directed=False)
        G.add_vertex(total_nodes)

        # Collect all edges
        edge_list = []

        # Step 1: Add base edges (1D line/ring)
        base_edges = self._create_1d_edges(
            n_nodes=self.base_nodes,
            start_vertex=0,
            periodic=self.periodic,
        )
        edge_list.extend(base_edges)

        # Step 2: Add fiber edges and anchor connections
        for base_idx in range(self.base_nodes):
            # Fiber starting vertex index
            fiber_start = self.base_nodes + base_idx * self.fiber_nodes

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
            'base_nodes': self.base_nodes,
            'fiber_nodes': self.fiber_nodes,
            'total_nodes': total_nodes,
            'structure': 'dirac_comb',
            'base_type': self.base_type,
            'periodic': self.periodic,
        }

        return G

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        # Total = base nodes + (base_nodes * fiber_nodes)
        return self.base_nodes * (1 + self.fiber_nodes)

    def get_base_vertex_indices(self) -> list:
        """Return list of base vertex indices."""
        return list(range(self.base_nodes))

    def get_fiber_vertex_indices(self, base_idx: int) -> list:
        """Return list of fiber vertex indices for a given base node.

        Parameters
        ----------
        base_idx : int
            Index of the base node

        Returns
        -------
        list
            List of vertex indices in the fiber attached to this base node
        """
        if base_idx < 0 or base_idx >= self.base_nodes:
            raise ValueError(f"base_idx must be in [0, {self.base_nodes})")
        fiber_start = self.base_nodes + base_idx * self.fiber_nodes
        return list(range(fiber_start, fiber_start + self.fiber_nodes))

    def __repr__(self) -> str:
        return (
            f"DiracCombGraphGT(base={self.base_nodes}, "
            f"fiber={self.fiber_nodes}, N={self.N})"
        )


# Backward compatibility alias
DiracCombGraph = DiracCombGraphGT
