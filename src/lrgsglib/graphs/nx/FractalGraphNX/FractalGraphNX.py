"""
FractalGraphNX: Abstract base class for fractal/self-similar graphs.

This module provides a unified interface for fractal graph models,
with common functionality for self-similar network structures.

Example
-------
>>> from lrgsglib.graphs.nx.DGMgraphNX import DGMgraphNX
>>> dgm = DGMgraphNX(n=5, pflip=0.1, seed=42)
>>> dgm.flip_random_fract_edges()
"""

from abc import abstractmethod
from typing import Any, Sequence
import random

import networkx as nx

from ....config.const import SG_GRAPHINT_REPR, COUNT_XERR_PATTERNS
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


class FractalGraphNX(SignedGraphNX):
    """
    Abstract base class for fractal/self-similar graph models.

    This class provides the common interface and functionality shared by
    all fractal graph types (DGM, Sierpinski, etc.).

    Parameters
    ----------
    n : int
        Iteration level or generation number for the fractal construction.
    sgpathn : str
        Subpath for storing graph data.
    std_fname : str
        Filename prefix for exports.
    with_positions : bool, default False
        Whether to compute and store node positions.
    only_const_mode : bool, default False
        If True, populate metadata only without building the graph.
    **kwargs
        Additional arguments passed to SignedGraphNX (pflip, seed, etc.).

    Attributes
    ----------
    n : int
        Iteration/generation level.
    syshape : int
        Number of nodes (depends on fractal formula).
    syshapePth : str
        Path string representation.
    H : nx.Graph
        Coordinate-labeled graph (if applicable).
    G : nx.Graph
        Integer-labeled graph.

    Notes
    -----
    Fractal graphs have node counts that grow according to specific formulas
    based on the iteration level n. Each fractal type has its own growth rate.
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
        self.n = n
        self.only_const_mode = only_const_mode
        self.with_positions = with_positions
        self.sgpathn = sgpathn
        self.std_fname = std_fname

        if not only_const_mode:
            self._init_network()
        else:
            self.syshape = self._compute_expected_nodes()
            self.syshapePth = f"N={self.syshape}"
            self.G = nx.Graph()

        super().__init__(self.G, **kwargs)

    @abstractmethod
    def _init_network(self) -> None:
        """Build the fractal graph at iteration level n."""
        raise NotImplementedError

    @abstractmethod
    def _compute_expected_nodes(self) -> int:
        """
        Compute expected number of nodes for iteration level n.

        Returns
        -------
        int
            Expected number of nodes.
        """
        raise NotImplementedError

    def get_expected_num_nodes(self) -> int:
        """
        Return the expected number of nodes for this fractal graph.

        Returns
        -------
        int
            Number of nodes.
        """
        return int(self.syshape)


class FractalNwContainerBase(dict):
    """
    Base class for fractal graph network container patterns.

    Parameters
    ----------
    graph : FractalGraphNX
        The parent fractal graph object.
    iterable : list, optional
        Initial keys to populate.
    constant : Any, optional
        Constant value to assign to initial keys.
    **kwargs
        Additional dict initialization arguments.
    """

    def __init__(
        self,
        graph: FractalGraphNX,
        iterable: Sequence = (),
        constant: Any = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.update((key, constant) for key in iterable)
        self.sg = graph
        self.rd = self.sg.graph_reprs
        self.rNodeFlip = {
            g: random.sample(list(self.sg.nodes_in(g)), self.sg.nflip)
            for g in self.rd
        }

        # Initialize common patterns
        self['rand'] = {g: list(self.sg.fleset[g]) for g in self.rd}
        self['randXERR'] = {
            g: self._get_rand_xerr_pattern(g)
            for g in self.rd
        }

    def get_links_XERR(self, node: Any, on_g: str = SG_GRAPHINT_REPR) -> list[tuple]:
        """Get star pattern: all edges incident to a node."""
        return [(node, nn) for nn in self.sg.get_graph_neighbors(node, on_g)]

    def _get_rand_xerr_pattern(self, on_g: str) -> list[tuple]:
        """Get random XERR pattern avoiding fully negative neighborhoods."""
        if COUNT_XERR_PATTERNS:
            return [
                k for i in self.rNodeFlip[on_g]
                for k in self.get_links_XERR(i, on_g)
            ]

        tmplst = list(self.rNodeFlip[on_g])
        grph = self.sg.gr[on_g]
        idx = 0
        patternList = []

        while idx < len(tmplst):
            leval = [
                all(nnn['weight'] == -1 for nnn in grph[nn].values())
                for nn in grph.neighbors(tmplst[idx])
            ]
            if any(leval):
                tmplst.pop(idx)
            else:
                patternList.extend(self.get_links_XERR(tmplst[idx], on_g))
                idx += 1

        return list(set(patternList))
