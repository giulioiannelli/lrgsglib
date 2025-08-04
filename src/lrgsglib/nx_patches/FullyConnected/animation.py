from typing import Union
import networkx as nx
import numpy as np
from ...plotlib import Colormap
from ..SignedGraph import SignedGraph
from ..funcs import signed_spectral_layout


def make_animation(
    fc: SignedGraph,
    fig,
    ax,
    frames,
    cmap: Union[str, Colormap] = "viridis",
):
    """Create animation for the fully connected graph."""
    G_nodecol = frames[0]
    G_edgecol = [
        "b" if (e[2]["weight"] > 0) else "r" for e in fc.G.edges(data=True)
    ]
    if fc.animation_graph_embedding == "sle":
        pos = signed_spectral_layout(fc.G)
    elif fc.animation_graph_embedding == "circular":
        pos = nx.circular_layout(fc.G)
    nodes = nx.draw_networkx_nodes(
        fc.G, pos=pos, node_color=G_nodecol, cmap=cmap
    )
    nx.draw_networkx_edges(fc.G, pos=pos, edge_color=G_edgecol)
    cbar = fig.colorbar(nodes)

    def animate(i):
        G_nodecol = frames[i]
        vmax = np.max(G_nodecol)
        vmin = np.min(G_nodecol)
        nx.draw_networkx_nodes(fc.G, pos=pos, node_color=G_nodecol, cmap=cmap)
        cbar.mappable.set_clim(vmin, vmax)

    return animate
