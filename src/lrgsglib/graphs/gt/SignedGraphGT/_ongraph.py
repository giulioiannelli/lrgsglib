"""Graph operations mixin for SignedGraphGT (NX-compatible API).

Methods here mirror ``graphs/nx/SignedGraphNX/_ongraph.py`` for feature
parity.  Methods already implemented directly on SignedGraphGT (flip_*,
unflip_all, get_random_links) are NOT duplicated here.
"""

from __future__ import annotations

import logging
from typing import Any, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ....config.const import SG_REPR, SG_ERRMSG_NFLIP
from ....config.errwar import NflipError

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraphGT import SignedGraphGT


def check_Ne_flips(self: "SignedGraphGT") -> None:
    """Validate that pflip produces at least one edge to flip."""
    ne_flips = int(self.pflip * self.num_edges)
    if ne_flips < 1:
        raise NflipError(SG_ERRMSG_NFLIP)


def set_edges_random_normal(
    self: "SignedGraphGT",
    mu: float = 1.0,
    sigma: float = 1.0,
    on_g: Optional[str] = None,
) -> None:
    """Set all edge weights from a normal distribution N(mu, sigma).

    Parameters
    ----------
    mu : float
        Mean of the distribution.
    sigma : float
        Standard deviation.
    on_g : str, optional
        Ignored (GT single representation).
    """
    sign_prop = self.G.edge_properties["sign"]
    for e in self.G.edges():
        sign_prop[e] = self._rng.normal(mu, sigma)
    self.invalidate_cache()


def load_vec_on_nodes(
    self: "SignedGraphGT",
    vec: NDArray,
    attr: str,
    on_g: Optional[str] = None,
) -> None:
    """Load a numpy vector as a vertex property map.

    Parameters
    ----------
    vec : ndarray, shape (N,)
        Values to assign to each vertex.
    attr : str
        Property map name (e.g. ``'eigV0'``).
    on_g : str, optional
        Ignored (GT single representation).
    """
    prop = self.G.new_vertex_property("double")
    for i, v in enumerate(self.G.vertices()):
        prop[v] = vec[i]
    self.G.vertex_properties[attr] = prop


def set_node_attributes(
    self: "SignedGraphGT",
    values: Any,
    attribute_name: str,
    on_g: Optional[str] = None,
) -> None:
    """Set vertex attributes from a dict or scalar.

    Parameters
    ----------
    values : dict or scalar
        If dict: ``{node_index: value}``.
        If scalar: applied to all vertices.
    attribute_name : str
        Property map name.
    on_g : str, optional
        Ignored (GT single representation).
    """
    if isinstance(values, dict):
        prop = self.G.new_vertex_property("double")
        for node, val in values.items():
            prop[self.G.vertex(node)] = val
        self.G.vertex_properties[attribute_name] = prop
    else:
        prop = self.G.new_vertex_property("double", val=float(values))
        self.G.vertex_properties[attribute_name] = prop


def get_random_edges_from_set(
    self: "SignedGraphGT",
    n: int = 1,
    edge_set: Optional[set] = None,
    on_g: Optional[str] = None,
) -> list:
    """Get n random edges from a specified set.

    Parameters
    ----------
    n : int
        Number of edges to sample.
    edge_set : set, optional
        Set of ``(u, v)`` tuples. If None, uses all edges.
    on_g : str, optional
        Ignored (GT single representation).

    Returns
    -------
    list of tuples
        Sampled edges.
    """
    if edge_set is None:
        edges = [(int(e.source()), int(e.target())) for e in self.G.edges()]
    else:
        edges = list(edge_set)

    if not edges:
        return []
    n = min(n, len(edges))
    indices = self._rng.choice(len(edges), size=n, replace=False)
    return [edges[i] for i in indices]
