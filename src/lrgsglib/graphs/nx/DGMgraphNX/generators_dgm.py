from networkx import Graph
import networkx as nx
from typing import Any

__all__ = ["dorogovtsev_goltsev_mendes_graph_FastPatch"]

def dorogovtsev_goltsev_mendes_graph_FastPatch(
    n: int,
    create_using: Any = None,
) -> Graph:
    """Generate a Dorogovtsev-Goltsev-Mendes graph.

    Parameters
    ----------
    n : int
        Iteration step for the construction. The resulting graph has
        ``3 * 2**n - 2`` nodes.
    create_using : Any, optional
        NetworkX graph type to create.

    Returns
    -------
    Graph
        A NetworkX graph instance.
    """
    return nx.dorogovtsev_goltsev_mendes_graph(n, create_using=create_using)
