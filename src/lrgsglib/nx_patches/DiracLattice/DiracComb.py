"""Dirac comb graph: 1D base with 1D fibers."""

from __future__ import annotations

import networkx as nx

from ..common import *
from .DiracLattice import DiracLatticeGraph
from ..MultispectralGraph.generators_msg import dirac_comb_graph


class DiracCombGraph(DiracLatticeGraph):
    """Dirac comb: 1D base network with 1D fiber networks.

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
    """

    def __init__(
        self,
        base_nodes: int,
        fiber_nodes: int,
        base_type: str = "line",
        periodic: bool = True,
        stdFnameSFFX: str = DCOMB_STDFN,
        sgpathn: str = DCOMB_SGPATH,
        **kwargs,
    ):
        self.base_nodes = base_nodes
        self.fiber_nodes = fiber_nodes
        self.base_type = base_type

        super().__init__(
            periodic=periodic,
            stdFnameSFFX=stdFnameSFFX,
            sgpathn=sgpathn,
            **kwargs
        )

    def _generate(self) -> nx.Graph:
        """Generate Dirac comb graph."""
        H, metadata = dirac_comb_graph(
            self.base_nodes,
            self.fiber_nodes,
            base_type=self.base_type,
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
        return self.base_nodes * (1 + self.fiber_nodes)

    def __repr__(self) -> str:
        return (
            f"DiracCombGraph(base={self.base_nodes}, "
            f"fiber={self.fiber_nodes}, N={self.N})"
        )
