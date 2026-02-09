"""
CompleteGraphGT: Abstract base class for complete graph models in graph-tool.

This mirrors the NX pattern where CompleteGraphNX is abstract and
FullyConnectedNX is the concrete implementation.
"""
from __future__ import annotations

from abc import abstractmethod
from typing import Optional

import numpy as np
import graph_tool.all as gt

from ..SignedGraphGT import SignedGraphGT


class CompleteGraphGT(SignedGraphGT):
    """
    Abstract base class for complete graph models using graph-tool.

    This class provides the common interface and functionality shared by
    all complete graph types (FullyConnected, Hopfield networks, etc.).

    For a concrete implementation, use FullyConnectedGT.

    Parameters
    ----------
    n : int
        Number of nodes in the complete graph.
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.
    with_positions : bool, default False
        Whether to compute and store node positions.
    **kwargs
        Additional arguments passed to SignedGraphGT.

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

    See Also
    --------
    FullyConnectedGT : Concrete implementation of complete graph.
    """

    def __init__(
        self,
        n: int,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        with_positions: bool = False,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self._N = n
        self.n = n
        self.only_const_mode = only_const_mode
        self.with_positions = with_positions
        self.syshape = n
        self.syshapePth = f"N={n}"

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        if not only_const_mode:
            G = self._init_network()
        else:
            G = gt.Graph(directed=False)

        super().__init__(G=G, pflip=pflip, seed=seed)

    @abstractmethod
    def _init_network(self) -> gt.Graph:
        """Build the complete graph. Must be implemented by subclasses."""
        raise NotImplementedError

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for this complete graph."""
        return self._N

    def get_expected_num_edges(self) -> int:
        """Return the expected number of edges for this complete graph."""
        return self._N * (self._N - 1) // 2


__all__ = ["CompleteGraphGT"]
