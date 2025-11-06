import logging
import random
import networkx as nx
from typing import Any, TYPE_CHECKING
from numpy.typing import NDArray

from ...config.const import SG_REPR, SG_ERRMSG_NFLIP
from ...config.errwar import NflipError

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def check_Ne_flips(self: "SignedGraph"):
    if self.Ne_flips < 1:
        raise NflipError(SG_ERRMSG_NFLIP)


def flip_sel_edges(
        self: "SignedGraph",
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


def flip_random_fract_edges(
        self: "SignedGraph",
        pflip: float = 0.0, 
        on_g: str = SG_REPR
) -> None:
    try:
        if pflip:
            self.pflip = pflip
            self.Ne_flips = int(self.pflip * self.Ne)
            self.check_Ne_flips()
            self.flip_sel_edges(
                self.get_random_links(self.Ne_flips, on_g), 
                on_g=on_g
            )
        else:
            self.Ne_flips = int(self.pflip * self.Ne)
            self.check_Ne_flips()
            self.flip_sel_edges(self.fleset[on_g], on_g=on_g)
    except NflipError:
        logger.error(SG_ERRMSG_NFLIP)
        pass


def unflip_all(self: "SignedGraph", on_g: str = SG_REPR):
    self.flip_sel_edges(1, on_g=on_g)


def set_edges_random_normal(
        self: "SignedGraph",
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
        self: "SignedGraph",
        vec: NDArray, 
        attr: str,
        on_g: str = SG_REPR
):
    vecNodeAttr = {nd: v for v, nd in zip(vec, self.gr[on_g].nodes)}
    nx.set_node_attributes(self.gr[on_g], values=vecNodeAttr, name=attr)


def load_eigV_on_graph(
        self: "SignedGraph",
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
        self: "SignedGraph",
        values: Any, 
        attribute_name: Any, 
        on_g: str = SG_REPR
):
    nx.set_node_attributes(self.gr[on_g], values, attribute_name)
