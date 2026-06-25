import logging
import random
import networkx as nx
import numpy as np
from typing import Any, Optional, TYPE_CHECKING
from numpy.typing import NDArray

from ....config.const import SG_REPR, SG_ERRMSG_NFLIP
from ....config.errwar import NflipError

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraphNX import SignedGraphNX


def check_Ne_flips(self: "SignedGraphNX"):
    if self.Ne_flips < 1:
        raise NflipError(SG_ERRMSG_NFLIP)


def flip_sel_edges(
        self: "SignedGraphNX",
        links: Any, 
        name: str = 'weight', 
        on_g: str = SG_REPR
) -> None:
    neg_weights_dict = {
        (u, v): -1 * self.get_edge_data(u, v, on_g=on_g) 
        for u, v in links
    }
    nx.set_edge_attributes(
        self.gr[on_g],
        values=neg_weights_dict, 
        name=name
    )

    self.fleset[on_g].update(links)
    self.lfeset[on_g].difference_update(links)

    self.upd_GraphRepr_All(on_g)
    self.upd_graph_matrices(on_g)


def get_random_edges_from_set(
        self: "SignedGraphNX",
        n: int = 1,
        edge_set: Optional[set] = None,
        on_g: str = SG_REPR
) -> list:
    """
    Get n random edges from a specified edge set efficiently.
    
    This method uses NumPy's random choice for efficient sampling from large
    edge sets without replacement. For large graphs, this is significantly
    faster than converting to list and using random.sample().
    
    Performance:
    - O(n) for conversion to array
    - O(k) for random selection where k is the sample size
    - Much faster than random.sample() for large edge sets
    
    Parameters
    ----------
    n : int, default 1
        Number of random edges to return.
    edge_set : set, optional
        The set of edges to sample from. If None, uses all edges (self.eset).
    on_g : str, default SG_REPR
        Graph representation to use.
        
    Returns
    -------
    list
        List of randomly selected edges as tuples.
    """
    if edge_set is None:
        edge_set = self.eset[on_g]
    
    # Convert set to numpy array for efficient indexing
    edge_array = np.array(list(edge_set), dtype=object)
    
    # Use numpy's random choice for efficient sampling without replacement
    # This is faster than random.sample for large sets
    indices = np.random.choice(len(edge_array), size=n, replace=False)
    
    # Return as list of tuples
    return [tuple(edge_array[i]) for i in indices]


def flip_random_fract_edges(
        self: "SignedGraphNX",
        pflip: Optional[float] = None, 
        on_g: str = SG_REPR
) -> None:
    try:
        if pflip:
            self.pflip = pflip
            self.Ne_flips = int(self.pflip * self.Ne)
            self.check_Ne_flips()
            # Use the new method to get random edges
            edges_to_flip = self.get_random_edges_from_set(
                self.Ne_flips, 
                self.eset[on_g], 
                on_g
            )
            self.flip_sel_edges(edges_to_flip, on_g=on_g)
        else:
            self.Ne_flips = int(self.pflip * self.Ne)
            self.check_Ne_flips()
            # Idempotent realize: SET the current flip set to negative weight
            # (-|w|) rather than toggling. A redundant call after the
            # construction-time auto-flip (disorder default) is then a no-op,
            # and the deferred (disorder=None) workflow converges here too.
            graph = self.gr[on_g]
            neg = {
                (u, v): -abs(self.get_edge_data(u, v, on_g=on_g))
                for (u, v) in self.fleset[on_g]
            }
            nx.set_edge_attributes(graph, values=neg, name='weight')
            self.upd_edge_sets(on_g)
            self.upd_GraphRepr_All(on_g)
            self.upd_graph_matrices(on_g)
    except NflipError:
        logger.error(SG_ERRMSG_NFLIP)
        pass


def unflip_all(self: "SignedGraphNX", on_g: str = SG_REPR):
    """
    Flip all negative edges back to positive.

    This reverses all edge sign flips, restoring the graph to all positive edges.

    Parameters
    ----------
    on_g : str, default SG_REPR
        Graph representation identifier.

    Notes
    -----
    This method flips all edges currently in fleset[on_g] (negative edges)
    back to positive, effectively restoring the base graph state.
    """
    # Create a copy of fleset to avoid modification during iteration
    edges_to_flip = set(self.fleset[on_g])
    if edges_to_flip:
        # Set all edges to positive weight (+1)
        import networkx as nx
        pos_weights_dict = {
            (u, v): 1.0
            for u, v in edges_to_flip
        }
        nx.set_edge_attributes(
            self.gr[on_g],
            values=pos_weights_dict,
            name='weight'
        )

        # Update the edge sets (opposite of flip_sel_edges)
        self.fleset[on_g].difference_update(edges_to_flip)
        self.lfeset[on_g].update(edges_to_flip)

        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices(on_g)


def set_edges_random_normal(
        self: "SignedGraphNX",
        mu: float = 1.0, 
        sigma: float = 1.0, 
        on_g: str = SG_REPR
):
    weights = {
        edge: random.normalvariate(mu, sigma) 
        for edge in self.gr[on_g].edges()
    }
    nx.set_edge_attributes(
        self.gr[on_g], 
        values=weights, 
        name='weight'
    )


def load_vec_on_nodes(
        self: "SignedGraphNX",
        vec: NDArray, 
        attr: str,
        on_g: str = SG_REPR
):
    vecNodeAttr = {nd: v for v, nd in zip(vec, self.gr[on_g].nodes)}
    nx.set_node_attributes(self.gr[on_g], values=vecNodeAttr, name=attr)


def load_eigV_on_graph(
        self: "SignedGraphNX",
        which: int = 0,
        on_g: str = SG_REPR,
        binarize: bool = False
):
    from . import _spectral

    if binarize:
        eigV = _spectral.get_eigV_bin_check(self, which=which)
    else:
        try:
            eigV = getattr(self, 'eigV', None)
            if eigV is None:
                raise AttributeError
            eigV = eigV[which]
        except (IndexError, AttributeError):
            _spectral.compute_k_eigvV(self, k=which + 1)
            eigV = getattr(self, 'eigV', None)
            if eigV is not None:
                eigV = eigV[which]

    if eigV is not None:
        load_vec_on_nodes(self, eigV, f"eigV{which}", on_g)


def set_node_attributes(
        self: "SignedGraphNX",
        values: Any, 
        attribute_name: Any, 
        on_g: str = SG_REPR
):
    nx.set_node_attributes(self.gr[on_g], values, attribute_name)
