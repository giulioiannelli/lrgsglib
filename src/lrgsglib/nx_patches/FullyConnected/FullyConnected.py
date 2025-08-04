from ...plotlib import Colormap
from ..common import *
from ..funcs import *
from ..SignedGraph.SignedGraph import SignedGraph
from .animation import make_animation as fc_make_animation
from typing import Union, Callable


class SignedGraph(SignedGraph):
    #
    def __init__(
        self,
        N: int = FC_N,
        with_positions: bool = False,
        mode_positions: Union[str, Callable] = "circular",
        anigemb: str = "sle",
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self._N = N
        self.with_positions = with_positions
        self.mode_positions = mode_positions
        self.only_const_mode = only_const_mode
        if not only_const_mode:
            self.__init_fullyconnected__()
        else:
            self.syshape = N
            self.syshapePth = f"N={N}"
            self.G = nx.Graph()
        super(SignedGraph, self).__init__(self.G, **kwargs)
        self.animation_graph_embedding = anigemb

    #
    def __init_fullyconnected__(self) -> None:
        self.G = nx.complete_graph(self._N)
        if self.with_positions:
            match self.mode_positions:
                case "circular":
                    pos = nx.circular_layout(self.G)
                case "sspectral":
                    pos = signed_spectral_layout(self.G)
                case callable(mode) if callable(mode):
                    pos = mode(self.G)
                case _:
                    raise ValueError(
                        f"Unsupported mode_positions: {self.mode_positions}"
                    )
            nx.set_node_attributes(self.G, pos, "pos")

    #
    def __repr__(self) -> str:
        return f"FullyConnected(N={self._N})"

    #
    def __str__(self) -> str:
        return f"FullyConnected(N={self._N})"

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the fully connected graph."""
        return self._N

    #
    def compute_hopfield_edges(self, mode: str) -> None:
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

    #
    def set_hopfield_edges(
        self, on_g: str = SG_GRAPH_REPR, **kwcompute_hopf_edges
    ) -> None:
        if not hasattr(self, "hopfield_edges"):
            self.compute_hopfield_edges(**kwcompute_hopf_edges)
        nx.set_edge_attributes(self.gr[on_g], self.hopfield_edges, "weight")

    #
    def make_animation(
        self, fig, ax, frames, cmap: Union[str, Colormap] = "viridis"
    ):
        return fc_make_animation(self, fig, ax, frames, cmap=cmap)
