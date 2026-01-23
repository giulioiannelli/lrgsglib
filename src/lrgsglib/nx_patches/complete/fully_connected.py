"""
FullyConnected: Backward-compatible fully connected graph wrapper.

This module provides the FullyConnected class that wraps the CompleteGraph
base class while maintaining full backward compatibility with the existing API.

Example
-------
>>> from lrgsglib.nx_patches.complete import FullyConnected
>>> fc = FullyConnected(N=50, pflip=0.1, seed=42)
>>> fc.flip_random_fract_edges()
>>> fc.N
50
"""

from typing import Callable, Union

import networkx as nx
import numpy as np

from ...config.const import FC_N, FC_PHTABB, FC_SGPATH, FC_STDFN, SG_REPR
from ...plotlib import Colormap
from ..funcs import signed_spectral_layout
from .complete_graph import CompleteGraph


class FullyConnected(CompleteGraph):
    """
    Signed fully connected (complete) graph.

    Parameters
    ----------
    N : int, default FC_N
        Number of nodes in the complete graph.
    sgpathn : str, default FC_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default FC_STDFN
        Filename suffix used in on-disk exports.
    with_positions : bool, default False
        Whether to store node positions for plotting.
    mode_positions : str | Callable, default "circular"
        Layout mode for positions. Options: "circular", "sspectral",
        or a callable that takes a graph and returns positions.
    anigemb : str, default "sle"
        Animation graph embedding mode.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraph`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Notes
    -----
    ``pflip`` defines how many edges are marked for sign flips, but weights
    are not set negative until you call ``flip_random_fract_edges`` or
    ``flip_sel_edges``.

    Examples
    --------
    >>> from lrgsglib.nx_patches.complete import FullyConnected
    >>> fc = FullyConnected(N=100, pflip=0.2, seed=1)
    >>> fc.flip_random_fract_edges()
    """

    def __init__(
        self,
        N: int = FC_N,
        sgpathn: str = FC_SGPATH,
        stdFnameSFFX: str = FC_STDFN,
        with_positions: bool = False,
        mode_positions: Union[str, Callable] = "circular",
        anigemb: str = "sle",
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.mode_positions = mode_positions
        self.animation_graph_embedding = anigemb
        self._init_std_fname(stdFnameSFFX)
        sgpath = FC_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=N,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            with_positions=with_positions,
            only_const_mode=only_const_mode,
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = FC_PHTABB + suffix

    def _init_network(self) -> None:
        """Build the fully connected graph."""
        self.G = nx.complete_graph(self._N)
        if self.with_positions:
            self._set_positions()

    def _set_positions(self) -> None:
        """Set node positions based on mode_positions."""
        match self.mode_positions:
            case "circular":
                pos = nx.circular_layout(self.G)
            case "sspectral":
                pos = signed_spectral_layout(self.G)
            case mode if callable(mode):
                pos = mode(self.G)
            case _:
                raise ValueError(
                    f"Unsupported mode_positions: {self.mode_positions}"
                )
        nx.set_node_attributes(self.G, pos, "pos")

    def __repr__(self) -> str:
        return f"FullyConnected(N={self._N})"

    def __str__(self) -> str:
        return f"FullyConnected(N={self._N})"

    def compute_hopfield_edges(self, mode: str) -> None:
        """
        Compute Hopfield-style edge weights.

        Parameters
        ----------
        mode : str
            Weight assignment mode:
            - "random": Random +1/-1 weights
            - "all+": All weights +1
            - "all-": All weights -1
        """
        if mode == "random":
            self.hopfield_edges = {
                (u, v): np.random.choice([-1, 1]) for u, v in self.G.edges()
            }
        elif mode == "all+":
            self.hopfield_edges = {(u, v): 1 for u, v in self.G.edges()}
        elif mode == "all-":
            self.hopfield_edges = {(u, v): -1 for u, v in self.G.edges()}
        else:
            raise ValueError(f"Unsupported mode: {mode}")

    def set_hopfield_edges(
        self, on_g: str = SG_REPR, **kwcompute_hopf_edges
    ) -> None:
        """
        Set Hopfield-style edge weights on the graph.

        Parameters
        ----------
        on_g : str, default 'G'
            Graph representation to modify.
        **kwcompute_hopf_edges
            Arguments passed to compute_hopfield_edges if not yet computed.
        """
        if not hasattr(self, "hopfield_edges"):
            self.compute_hopfield_edges(**kwcompute_hopf_edges)
        nx.set_edge_attributes(self.gr[on_g], self.hopfield_edges, "weight")

    def make_animation(
        self, fig, ax, frames, cmap: Union[str, Colormap] = "viridis"
    ):
        """Create animation from frames data."""
        from ..FullyConnected.animation import make_animation as fc_make_animation
        return fc_make_animation(self, fig, ax, frames, cmap=cmap)
