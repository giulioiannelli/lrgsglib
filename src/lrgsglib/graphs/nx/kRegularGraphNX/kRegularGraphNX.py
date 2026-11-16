"""
kRegularGraphNX: Random k-regular graph.

A k-regular graph is a graph where every node has exactly degree k.
Useful as a baseline for studying disorder effects since the degree
distribution has zero variance.

``method='native'`` uses ``networkx.random_regular_graph``. ``method='pairing'``
uses the exact uniform pairing sampler shared with the GT engine, so the same
seed gives the same graph on both. The pairing draw takes its own child of the
seed (``pairing_stream``) and does not use the process-global RNGs, but the
``SignedGraphNX`` constructor reseeds them for either method.

Example
-------
>>> from lrgsglib.graphs.nx import kRegularGraphNX
>>> kreg = kRegularGraphNX(n=100, k=4, pflip=0.1, seed=42)
>>> kreg.flip_random_fract_edges()
"""

from typing import Optional

import networkx as nx

from ....config.const import (
    KREG_METHOD,
    KREG_METHOD_PAIRING,
    KREG_METHODS,
    KREG_PHTABB,
)
from ..._shared._regular_pairing import pairing_stream, random_regular_edges
from ..RandomGraphNX.RandomGraphNX import RandomGraphNX

# Constants
KREG_SGPATH = ""
KREG_STDFN = ""


class kRegularGraphNX(RandomGraphNX):
    """
    Signed random k-regular graph.

    Every node has exactly degree k. The topology comes from
    ``networkx.random_regular_graph`` (``method='native'``) or from the
    exact uniform pairing sampler (``method='pairing'``).

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Degree of each node. Must satisfy n*k even (for undirected graph).
    sgpathn : str, default KREG_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default KREG_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    method : {'native', 'pairing'}, default 'native'
        Topology sampler: networkx's generator, or the exact uniform pairing
        sampler (same graph as GT for the same seed).
    **kwargs : Any
        Forwarded to SignedGraphNX (pflip, seed, etc.). ``seed`` also fixes
        the topology: ``method='native'`` passes it to networkx,
        ``method='pairing'`` draws from ``pairing_stream(seed)``.

    Notes
    -----
    k-regular graphs are useful for:
    - Baseline comparisons (zero degree variance)
    - Studying pure structural effects without degree heterogeneity
    - Expander graph constructions

    Examples
    --------
    >>> from lrgsglib.graphs.nx import kRegularGraphNX
    >>> kreg = kRegularGraphNX(n=100, k=6, pflip=0.1, seed=42)
    >>> kreg.flip_random_fract_edges()
    >>> all(d == 6 for _, d in kreg.G.degree())
    True
    """

    def __init__(
        self,
        n: int,
        k: int,
        sgpathn: str = KREG_SGPATH,
        stdFnameSFFX: str = KREG_STDFN,
        only_const_mode: bool = False,
        method: str = KREG_METHOD,
        **kwargs,
    ) -> None:
        if (n * k) % 2 != 0:
            raise ValueError(f"n*k must be even, got n={n}, k={k}, n*k={n*k}")
        if k >= n:
            raise ValueError(f"k must be less than n, got k={k}, n={n}")
        if method not in KREG_METHODS:
            raise ValueError(f"method must be one of {KREG_METHODS}")

        self.k = k
        self.method = method
        self._init_std_fname(stdFnameSFFX)
        sgpath = KREG_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,  # k-regular graphs are typically connected for k >= 2
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        self.std_fname = KREG_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate random k-regular graph (seeded via ``_rng_seed``)."""
        seed = getattr(self, "_rng_seed", None)
        if self.method == KREG_METHOD_PAIRING:
            return self._graph_from_pairing(seed)
        return nx.random_regular_graph(self.k, self.n, seed=seed)

    def _graph_from_pairing(self, seed: Optional[int]) -> nx.Graph:
        """Undirected NX graph of the exact pairing sampler's edge list.

        Drawn from ``pairing_stream(seed)``, so the edges equal
        ``random_regular_edges(n, k, pairing_stream(seed))``.
        """
        edges = random_regular_edges(self.n, self.k, pairing_stream(seed))
        G = nx.Graph()  # noqa: raw-graph
        G.add_nodes_from(range(self.n))
        G.add_edges_from(edges.tolist())
        return G

    def _compute_syshapePth(self) -> str:
        path = f"N={self.n}_k={self.k}"
        if self.method != KREG_METHOD:
            path += f"_method={self.method}"
        return path
