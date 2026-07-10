"""
RandomGraphNX: Abstract base class for random graphs.

This module provides a unified interface for random graph models,
with common functionality shared across Erdos-Renyi, Watts-Strogatz,
Barabasi-Albert, and other random graph types.

Example
-------
>>> from lrgsglib.graphs.nx import ErdosRenyiNX
>>> er = ErdosRenyiNX(n=100, p=0.05, pflip=0.1, seed=42)
>>> er.flip_random_fract_edges()
>>> er.N
# Returns number of nodes in giant component
"""

from abc import abstractmethod
from typing import Any, Optional, Sequence
import random

import networkx as nx

from ....config.const import SG_GRAPHINT_REPR, COUNT_XERR_PATTERNS
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


class RandomGraphNX(SignedGraphNX):
    """
    Abstract base class for random graph models.

    This class provides the common interface and functionality shared by
    all random graph types (Erdos-Renyi, Watts-Strogatz, Barabasi-Albert, etc.).

    Parameters
    ----------
    n : int
        Number of nodes in the graph (or initial graph before giant component extraction).
    sgpathn : str
        Subpath for storing graph data.
    std_fname : str
        Filename prefix for exports.
    only_const_mode : bool, default False
        If True, populate metadata only without building the graph.
    extract_giant_component : bool, default True
        If True, extract and use only the largest connected component.
    **kwargs
        Additional arguments passed to SignedGraphNX (pflip, seed, etc.).

    Attributes
    ----------
    n : int
        Initial number of nodes requested.
    syshape : int
        Actual number of nodes (may differ after giant component extraction).
    syshapePth : str
        Path string representation of system parameters.

    Notes
    -----
    Subclasses must implement:
    - `_generate_graph()`: Generate the random graph using appropriate algorithm
    - `_compute_syshapePth()`: Return path string for this graph configuration
    """

    def __init__(
        self,
        n: int,
        sgpathn: str = "",
        std_fname: str = "",
        only_const_mode: bool = False,
        extract_giant_component: bool = True,
        **kwargs,
    ) -> None:
        self.n = n
        self.only_const_mode = only_const_mode
        self.extract_giant_component = extract_giant_component
        self.sgpathn = sgpathn
        self.std_fname = std_fname
        # Capture the seed BEFORE the graph is drawn: _generate_graph runs ahead
        # of SignedGraphNX.__init__, so subclasses reading a seed attribute
        # (e.g. ExtendedBarabasiAlbertNX's np.random.default_rng(self._rng_seed))
        # would otherwise always see None and draw irreproducibly.
        self._rng_seed = kwargs.get("seed")

        if not only_const_mode:
            self._init_network()
        else:
            self.G = nx.Graph()
            self.syshape = n
            self.syshapePth = self._compute_syshapePth()

        super().__init__(self.G, **kwargs)

    def _init_network(self) -> None:
        """Initialize the network graph."""
        G = self._generate_graph()

        if self.extract_giant_component and not nx.is_connected(G):
            CC = nx.connected_components(G)
            GC = max(CC, key=len)
            # Relabel the giant component to contiguous 0..N-1 integers (same
            # invariant as ErdosRenyiNX): dropping the non-giant nodes otherwise
            # leaves gaps in the labels, which breaks every positional consumer
            # (the statsys sweepers index the state array by node id). Node
            # attributes (e.g. 'pos') travel with the relabel.
            self.G = nx.convert_node_labels_to_integers(G.subgraph(GC).copy())
        else:
            self.G = G

        self.syshape = self.G.number_of_nodes()
        self.syshapePth = self._compute_syshapePth()

    @abstractmethod
    def _generate_graph(self) -> nx.Graph:
        """
        Generate the random graph.

        Returns
        -------
        nx.Graph
            The generated graph.

        Notes
        -----
        Subclasses must implement this method to generate their specific
        graph type using the appropriate NetworkX function or custom algorithm.
        """
        raise NotImplementedError

    @abstractmethod
    def _compute_syshapePth(self) -> str:
        """
        Compute the path string representation.

        Returns
        -------
        str
            Path string encoding graph parameters (e.g., "N=100_p=0.05").
        """
        raise NotImplementedError

    def get_expected_num_nodes(self) -> int:
        """
        Return the expected number of nodes for this random graph.

        Returns
        -------
        int
            Number of nodes in the graph.
        """
        return int(self.syshape)


class RandomNwContainerBase(dict):
    """
    Base class for random graph network container patterns.

    This provides common functionality for all random graph
    nwContainer classes, including star patterns and random edge selection.

    Parameters
    ----------
    graph : RandomGraphNX
        The parent random graph object.
    iterable : list, optional
        Initial keys to populate.
    constant : Any, optional
        Constant value to assign to initial keys.
    **kwargs
        Additional dict initialization arguments.
    """

    def __init__(
        self,
        graph: RandomGraphNX,
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
        """
        Get star pattern: all edges incident to a node.

        Parameters
        ----------
        node : Any
            The center node.
        on_g : str
            Graph representation ('G' or 'H').

        Returns
        -------
        list[tuple]
            List of edges as (node, neighbor) tuples.
        """
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
