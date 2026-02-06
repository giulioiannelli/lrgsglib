"""
CompleteGraph: Abstract base class for fully connected graphs.

This module provides a unified interface for complete graph models,
with common functionality for fully connected networks.

Example
-------
>>> from lrgsglib.nx_patches.complete import FullyConnected
>>> fc = FullyConnected(N=50, pflip=0.1, seed=42)
>>> fc.flip_random_fract_edges()
>>> fc.N
50
"""

from abc import abstractmethod
from typing import Any, Callable, Union

import networkx as nx

from ...config.const import SG_GRAPHINT_REPR
from ..SignedGraph.SignedGraph import SignedGraph


class CompleteGraph(SignedGraph):
    """
    Abstract base class for fully connected graph models.

    This class provides the common interface and functionality shared by
    all complete graph types (FullyConnected, Hopfield networks, etc.).

    Parameters
    ----------
    n : int
        Number of nodes in the complete graph.
    sgpathn : str
        Subpath for storing graph data.
    std_fname : str
        Filename prefix for exports.
    with_positions : bool, default False
        Whether to compute and store node positions.
    only_const_mode : bool, default False
        If True, populate metadata only without building the graph.
    **kwargs
        Additional arguments passed to SignedGraph (pflip, seed, etc.).

    Attributes
    ----------
    n : int
        Number of nodes.
    syshape : int
        System shape (equals n for complete graphs).
    syshapePth : str
        Path string representation.

    Notes
    -----
    Complete graphs have n*(n-1)/2 edges for n nodes.
    """

    def __init__(
        self,
        n: int,
        sgpathn: str = "",
        std_fname: str = "",
        with_positions: bool = False,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self._N = n
        self.only_const_mode = only_const_mode
        self.with_positions = with_positions
        self.sgpathn = sgpathn
        self.std_fname = std_fname
        self.syshape = n
        self.syshapePth = f"N={n}"

        if not only_const_mode:
            self._init_network()
        else:
            self.G = nx.Graph()

        super().__init__(self.G, **kwargs)

    @abstractmethod
    def _init_network(self) -> None:
        """Build the complete graph."""
        raise NotImplementedError

    def get_expected_num_nodes(self) -> int:
        """
        Return the expected number of nodes for this complete graph.

        Returns
        -------
        int
            Number of nodes.
        """
        return self._N

    def get_expected_num_edges(self) -> int:
        """
        Return the expected number of edges for this complete graph.

        Returns
        -------
        int
            Number of edges (n*(n-1)/2).
        """
        return self._N * (self._N - 1) // 2
