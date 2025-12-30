"""Multiplicative cascade graph generator."""

from __future__ import annotations

import networkx as nx

from ..common import *
from ..MultispectralGraph import MultispectralGraph
from ..MultispectralGraph.generators_msg import (
    multiplicative_cascade_graph,
    multiplicative_cascade_probability_matrix,
    multiplicative_cascade_exp_clocks,
)


class MultiplicativeCascadeGraph(MultispectralGraph):
    """Multiplicative cascade graph with hierarchical structure.

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
        exp_clocks is 2-16x faster depending on graph size
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
        stdFnameSFFX: str = MC_STDFN,
        sgpathn: str = MC_SGPATH,
        **kwargs,
    ):
        # Store parameters
        self.p1, self.p2, self.p3, self.p4 = p1, p2, p3, p4
        self.fraction = fraction
        self.iterations = iterations
        self.stochastic = stochastic
        self.periodic = periodic
        self.variant = variant

        # Computed attributes
        self.seed_probabilities = (p1, p2, p3, p4)
        self.sample_fraction = fraction
        self.probability_matrix = None

        super().__init__(
            stdFnameSFFX=stdFnameSFFX,
            sgpathn=sgpathn,
            **kwargs
        )

    def _generate(self) -> nx.Graph:
        """Generate multiplicative cascade graph."""
        # Compute probability matrix
        self.probability_matrix = multiplicative_cascade_probability_matrix(
            self.p1, self.p2, self.p3, self.p4,
            iterations=self.iterations,
            stochastic=self.stochastic,
        )

        # Generate graph using appropriate variant
        if self.variant == "exp_clocks":
            H = multiplicative_cascade_exp_clocks(
                self.probability_matrix,
                self.fraction,
                periodic=self.periodic
            )
        else:
            H = multiplicative_cascade_graph(
                self.p1, self.p2, self.p3, self.p4,
                fraction=self.fraction,
                iterations=self.iterations,
                stochastic=self.stochastic,
                probabilities=self.probability_matrix,
                periodic=self.periodic,
            )

        self.H = H
        return H

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        if self.only_const_mode or self.probability_matrix is None:
            return 0
        size = int(self.probability_matrix.shape[0])
        return max(1, int(round(self.sample_fraction * (size ** 2))))

    def __repr__(self) -> str:
        return f"MultiplicativeCascadeGraph(N={self.N}, variant={self.variant})"
