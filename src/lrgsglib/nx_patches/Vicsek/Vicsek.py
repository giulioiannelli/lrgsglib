"""Vicsek (Palla-Lovász-Vicsek) graph generator."""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np

from ..common import *
from ..MultispectralGraph import MultispectralGraph
from ..MultispectralGraph.generators_msg import (
    palla_lovasz_vicksek_graph,
    initial_measure,
)


class VicsekGraph(MultispectralGraph):
    """Vicsek (Palla-Lovász-Vicsek) graph via Kronecker product.

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
    """

    def __init__(
        self,
        N: int,
        k: int,
        pij: np.ndarray | None = None,
        m: int = 3,
        symmetric: bool = True,
        stdFnameSFFX: str = VCK_STDFN,
        sgpathn: str = VCK_SGPATH,
        **kwargs,
    ):
        # Store parameters
        self.N_target = N
        self.k = k
        self.pij_input = pij
        self.m = m
        self.symmetric = symmetric

        # Computed attributes
        self.probability_matrix = None
        self.sample_fraction = 1.0

        super().__init__(
            stdFnameSFFX=stdFnameSFFX,
            sgpathn=sgpathn,
            **kwargs
        )

    def _generate(self) -> nx.Graph:
        """Generate Vicsek graph."""
        # Generate or use provided initial measure
        if self.pij_input is None:
            pij = initial_measure(self.m, symmetric=self.symmetric)
        else:
            pij = self.pij_input

        # Generate graph
        H, probability_matrix = palla_lovasz_vicksek_graph(
            self.N_target, pij, self.k
        )

        self.probability_matrix = probability_matrix
        self.H = H
        return H

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        if self.only_const_mode:
            return 0
        return self.N_target

    def __repr__(self) -> str:
        return f"VicsekGraph(N={self.N}, k={self.k})"
