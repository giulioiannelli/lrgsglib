"""
Unified network container mixin for graph types.

Provides a standardized base class for the nwContainer pattern used
across different graph families for edge pattern selection.
"""
import random
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from ..SignedGraph.SignedGraph import SignedGraph

from ...config.const import SG_REPR, COUNT_XERR_PATTERNS


class nwContainerMixin(dict):
    """
    Base class for network container patterns.

    The nwContainer pattern provides standardized edge selection patterns
    for use in dynamics simulations (Ising, Contact Process, etc.). Each
    graph type can extend this base to add geometry-specific patterns.

    This mixin provides common functionality:
    - Random node selection based on pflip
    - XERR (star) pattern: all edges from a node
    - Base 'rand' pattern: all flip-marked edges

    Subclasses should add geometry-specific patterns like:
    - ZERR (loop) patterns for lattices
    - Cluster-based patterns for hierarchical graphs

    Parameters
    ----------
    sg : SignedGraph
        The signed graph instance this container belongs to.
    iterable : list, optional
        Keys to initialize with a constant value.
    constant : Any, optional
        Constant value for initialized keys.
    **kwargs
        Additional dict initialization arguments.

    Attributes
    ----------
    sg : SignedGraph
        Reference to the parent signed graph.
    rd : list
        List of graph representation keys.
    rNodeFlip : dict
        Random nodes selected for pattern generation.

    Examples
    --------
    >>> class Lattice2DContainer(nwContainerMixin):
    ...     def _init_patterns(self):
    ...         super()._init_patterns()
    ...         self['singleZERR'] = self._get_zerr_patterns()
    """

    sg: "SignedGraph"
    rd: list
    rNodeFlip: dict

    def __init__(
        self,
        sg: "SignedGraph",
        iterable: list = None,
        constant: Any = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        if iterable:
            self.update((key, constant) for key in iterable)

        self.sg = sg
        self.rd = self.sg.graph_reprs
        self._init_random_nodes()
        self._init_patterns()

    def _init_random_nodes(self) -> None:
        """Initialize random node selection based on nflip."""
        self.rNodeFlip = {
            g: random.sample(list(self.sg.gr[g].nodes()), self.sg.nflip)
            for g in self.rd
        }

    def _init_patterns(self) -> None:
        """
        Initialize standard patterns.

        Override this method in subclasses to add geometry-specific patterns.
        Always call super()._init_patterns() first to get base patterns.
        """
        # Base patterns available for all graph types
        self["rand"] = {g: list(self.sg.fleset[g]) for g in self.rd}
        self["randXERR"] = {g: self.get_rand_pattern("XERR", on_g=g) for g in self.rd}

    def get_links_XERR(self, node: Any, on_g: str = SG_REPR) -> list:
        """
        Get star pattern (all edges from a node).

        Parameters
        ----------
        node : Any
            The center node of the star.
        on_g : str, default SG_REPR
            Graph representation to use.

        Returns
        -------
        list
            List of (node, neighbor) tuples for all edges from the node.
        """
        return [(node, nn) for nn in self.sg.get_graph_neighbors(node, on_g)]

    def get_rand_pattern(self, mode: str, on_g: str = SG_REPR) -> list:
        """
        Get random pattern of specified type.

        Parameters
        ----------
        mode : str
            Pattern type: 'XERR' for star patterns.
        on_g : str, default SG_REPR
            Graph representation to use.

        Returns
        -------
        list
            List of edges forming the pattern.
        """
        match mode:
            case "XERR":
                if COUNT_XERR_PATTERNS:
                    pattern_list = [
                        k
                        for i in self.rNodeFlip[on_g]
                        for k in self.get_links_XERR(i, on_g)
                    ]
                else:
                    pattern_list = self._get_filtered_xerr_pattern(on_g)
            case _:
                pattern_list = []

        return list(set(pattern_list))

    def _get_filtered_xerr_pattern(self, on_g: str) -> list:
        """
        Get XERR pattern with filtering for adjacent negative edges.

        Parameters
        ----------
        on_g : str
            Graph representation to use.

        Returns
        -------
        list
            Filtered list of edges.
        """
        tmp_lst = list(self.rNodeFlip[on_g])
        grph = self.sg.gr[on_g]
        idx = 0
        pattern_list = []

        while idx < len(tmp_lst):
            node = tmp_lst[idx]
            # Check if any neighbor has all negative edges
            has_all_neg_neighbor = any(
                all(nnn.get("weight", 1) == -1 for nnn in grph[nn].values())
                for nn in grph.neighbors(node)
            )

            if has_all_neg_neighbor:
                tmp_lst.pop(idx)
            else:
                pattern_list.extend(self.get_links_XERR(node, on_g))
                idx += 1

        return pattern_list


class LatticenWContainerMixin(nwContainerMixin):
    """
    Extended nwContainer for lattice geometries.

    Adds ZERR (loop) patterns specific to different lattice geometries
    (triangular, squared, hexagonal).
    """

    def _init_patterns(self) -> None:
        """Initialize lattice-specific patterns."""
        super()._init_patterns()

        # Add central edge patterns
        if hasattr(self.sg, "get_central_edge"):
            self.centedge = {g: self.sg.get_central_edge(g) for g in self.rd}
            self["single"] = {g: [self.centedge[g]] for g in self.rd}
            self["singleXERR"] = {
                g: self.get_links_XERR(self.centedge[g][0], g) for g in self.rd
            }

    def get_links_rball(
        self,
        R: int = 1,
        center: Any = None,
        on_g: str = SG_REPR,
    ) -> set:
        """
        Get edges within radius R of a center node.

        Parameters
        ----------
        R : int, default 1
            Radius of the ball.
        center : Any, optional
            Center node. If None, uses central edge node.
        on_g : str, default SG_REPR
            Graph representation to use.

        Returns
        -------
        set
            Set of edges within the ball.
        """
        from ..funcs.neighbors import get_neighbors_at_distance

        graph = self.sg.gr[on_g]
        if center is None:
            center = self.centedge[on_g][0]

        neighs_to_flip = get_neighbors_at_distance(graph, center, R)
        links = {
            (node, neighbor)
            for node in neighs_to_flip
            for neighbor in graph.neighbors(node)
        }
        return links
