"""
kRegularGraph: Random k-regular graph.

A k-regular graph is a graph where every node has exactly degree k.
Useful as a baseline for studying disorder effects since the degree
distribution has zero variance.

Example
-------
>>> from lrgsglib.nx_patches.random import kRegularGraph
>>> kreg = kRegularGraph(n=100, k=4, pflip=0.1, seed=42)
>>> kreg.flip_random_fract_edges()
"""

import networkx as nx

from .random_graph import RandomGraph


# Constants
KREG_PHTABB = "kreg"
KREG_SGPATH = ""
KREG_STDFN = ""


class kRegularGraph(RandomGraph):
    """
    Signed random k-regular graph.

    Every node has exactly degree k. The graph is generated using
    the configuration model with a constant degree sequence.

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
    **kwargs : Any
        Forwarded to SignedGraph (pflip, seed, etc.).

    Notes
    -----
    k-regular graphs are useful for:
    - Baseline comparisons (zero degree variance)
    - Studying pure structural effects without degree heterogeneity
    - Expander graph constructions

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import kRegularGraph
    >>> kreg = kRegularGraph(n=100, k=6, pflip=0.1, seed=42)
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
        **kwargs,
    ) -> None:
        if (n * k) % 2 != 0:
            raise ValueError(f"n*k must be even, got n={n}, k={k}, n*k={n*k}")
        if k >= n:
            raise ValueError(f"k must be less than n, got k={k}, n={n}")

        self.k = k
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
        """Generate random k-regular graph."""
        return nx.random_regular_graph(self.k, self.n)

    def _compute_syshapePth(self) -> str:
        return f"N={self.n}_k={self.k}"
