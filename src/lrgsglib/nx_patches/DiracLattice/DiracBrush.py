"""Dirac brush graph: 2D base with 1D fibers."""

from __future__ import annotations

import networkx as nx

from ..common import *
from .DiracLattice import DiracLatticeGraph
from ..MultispectralGraph.generators_msg import dirac_brush_graph


class DiracBrushGraph(DiracLatticeGraph):
    """Dirac brush: 2D base network with 1D fiber networks.

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
    """

    def __init__(
        self,
        base_x: int,
        base_y: int,
        fiber_nodes: int,
        periodic: bool = True,
        stdFnameSFFX: str = DBRUSH_STDFN,
        sgpathn: str = DBRUSH_SGPATH,
        **kwargs,
    ):
        self.base_x = base_x
        self.base_y = base_y
        self.fiber_nodes = fiber_nodes

        super().__init__(
            periodic=periodic,
            stdFnameSFFX=stdFnameSFFX,
            sgpathn=sgpathn,
            **kwargs
        )

    def _generate(self) -> nx.Graph:
        """Generate Dirac brush graph."""
        H, metadata = dirac_brush_graph(
            self.base_x,
            self.base_y,
            self.fiber_nodes,
            periodic=self.periodic
        )

        # Store Dirac metadata
        self.dirac_structure = metadata
        self.base_graph = metadata['base_graph']
        self.fiber_graph = metadata['fiber_graph']
        self.H = H

        return H

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        if self.only_const_mode:
            return 0
        # Total = base nodes + (base_nodes × fiber_nodes)
        base_nodes = self.base_x * self.base_y
        return base_nodes * (1 + self.fiber_nodes)

    def __repr__(self) -> str:
        return (
            f"DiracBrushGraph(base={self.base_x}x{self.base_y}, "
            f"fiber={self.fiber_nodes}, N={self.N})"
        )
