"""
ErdosRenyi: Backward-compatible Erdos-Renyi random graph wrapper.

This module provides the ErdosRenyi class that wraps the RandomGraph base class
while maintaining full backward compatibility with the existing API.

Example
-------
>>> from lrgsglib.nx_patches.random import ErdosRenyi
>>> er = ErdosRenyi(n=200, p=0.05, pflip=0.2, seed=42)
>>> er.flip_random_fract_edges()
>>> er.N
# Returns number of nodes in giant component
"""

from typing import Any
import random

import networkx as nx

from ...config.const import (
    COUNT_XERR_PATTERNS,
    ER_ONREP,
    ER_PHTABB,
    ER_SGPATH,
    ER_STDFN,
)
from ..funcs import get_smallest_cycle_graph_node
from ..SignedGraph.SignedGraph import SignedGraph
from .random_graph import RandomGraph, RandomNwContainerBase


class ErdosRenyi(RandomGraph):
    """
    Signed Erdos-Renyi graph generator (largest connected component).

    Parameters
    ----------
    n : int
        Number of nodes in the initial Erdos-Renyi graph.
    p : float
        Edge probability for the Erdos-Renyi model.
    sgpathn : str, default ER_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default ER_STDFN
        Filename suffix used in on-disk exports.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraph`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Notes
    -----
    The graph is generated with NetworkX and then reduced to its largest
    connected component before initializing ``SignedGraph``. This ensures
    that most downstream algorithms operate on a connected graph.

    ``pflip`` defines how many edges are marked for sign flips, but weights
    are not set negative until you call ``flip_random_fract_edges`` or
    ``flip_sel_edges``, unless weights already exist on the input graph.

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import ErdosRenyi
    >>> er = ErdosRenyi(n=200, p=0.05, pflip=0.2, seed=1)
    >>> er.flip_random_fract_edges()  # apply sign flips
    """

    def __init__(
        self,
        n: int,
        p: float,
        sgpathn: str = ER_SGPATH,
        stdFnameSFFX: str = ER_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.p = p
        self._init_std_fname(stdFnameSFFX)
        sgpath = ER_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=True,
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = ER_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate Erdos-Renyi random graph."""
        return nx.erdos_renyi_graph(self.n, self.p)

    def _compute_syshapePth(self) -> str:
        """Compute path string for Erdos-Renyi parameters."""
        return f"N={self.n}_p={self.p:.3g}"

    class nwContainer(RandomNwContainerBase):
        """Network container for Erdos-Renyi patterns."""

        def __init__(
            self,
            er: SignedGraph,
            iterable: list = None,
            constant: Any = None,
            **kwargs,
        ):
            if iterable is None:
                iterable = []
            # Use 'er' as the graph reference for backward compatibility
            self.er = er
            super().__init__(er, iterable, constant, **kwargs)

            # Add ZERR patterns
            self['randZERR'] = {
                g: self.get_rand_pattern('ZERR', on_g=g)
                for g in self.rd
            }

        def get_links_ZERR(self, node: Any, on_g: str = ER_ONREP) -> list:
            """Get smallest cycle pattern for zero-error."""
            smallest_cycle = get_smallest_cycle_graph_node(self.sg.gr[on_g], node)
            edges = []
            if smallest_cycle:
                edges = [
                    (smallest_cycle[i], smallest_cycle[i + 1])
                    for i in range(len(smallest_cycle) - 1)
                ]
            return edges

        def get_rand_pattern(self, mode: str, on_g: str = ER_ONREP) -> list:
            """Get random pattern of given mode."""
            if mode == "XERR":
                return self._get_rand_xerr_pattern(on_g)
            elif mode == "ZERR":
                return self._get_rand_zerr_pattern(on_g)
            return []

        def _get_rand_zerr_pattern(self, on_g: str) -> list:
            """Get random ZERR pattern based on smallest cycles."""
            return [
                links
                for node in self.rNodeFlip[on_g]
                if (links := self.get_links_ZERR(node, on_g))
            ]
