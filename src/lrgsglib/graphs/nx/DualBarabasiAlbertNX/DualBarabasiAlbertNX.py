"""
DualBarabasiAlbertNX: Dual BA model with two attachment modes.

NetworkX dual BA model where new nodes attach with m1 edges (probability p)
or m2 edges (probability 1-p) via preferential attachment.

Example
-------
>>> from lrgsglib.graphs.nx import DualBarabasiAlbertNX
>>> g = DualBarabasiAlbertNX(n=500, m1=2, m2=4, p=0.5, seed=42)
"""

import networkx as nx

from ..RandomGraphNX.RandomGraphNX import RandomGraphNX


# Constants
PA_PHTABB = "pa"
PA_SGPATH = ""
PA_STDFN = ""


class DualBarabasiAlbertNX(RandomGraphNX):
    """
    Dual Barabasi-Albert model with two attachment modes.

    NetworkX dual BA model where new nodes attach with m1 edges (probability p)
    or m2 edges (probability 1-p) via preferential attachment.

    Parameters
    ----------
    n : int
        Number of nodes.
    m1 : int
        Number of edges in mode 1 (with probability p).
    m2 : int
        Number of edges in mode 2 (with probability 1-p).
    p : float, default 0.5
        Probability of using m1 edges (vs m2 edges).
    sgpathn : str, default PA_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default PA_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraphNX (pflip, seed, etc.).

    Notes
    -----
    - p=1: All nodes use m1 edges
    - p=0: All nodes use m2 edges
    - 0<p<1: Mixed (heterogeneous attachment)

    This models networks where different nodes have different "budgets"
    for forming connections.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import DualBarabasiAlbertNX
    >>> # Half use 2 edges, half use 4 edges
    >>> g = DualBarabasiAlbertNX(n=500, m1=2, m2=4, p=0.5, seed=42)
    """

    def __init__(
        self,
        n: int,
        m1: int,
        m2: int,
        p: float = 0.5,
        sgpathn: str = PA_SGPATH,
        stdFnameSFFX: str = PA_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if m1 < 1 or m2 < 1:
            raise ValueError(f"m1 and m2 must be >= 1, got m1={m1}, m2={m2}")
        if not 0 <= p <= 1:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if max(m1, m2) >= n:
            raise ValueError(f"max(m1, m2) must be < n")

        self.m1 = m1
        self.m2 = m2
        self.p = p
        self._init_std_fname(stdFnameSFFX)
        sgpath = PA_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=False,
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = PA_PHTABB + "_dual" + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate dual BA graph."""
        seed = getattr(self, '_rng_seed', None)
        return nx.dual_barabasi_albert_graph(self.n, self.m1, self.m2, self.p, seed=seed)

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"N={self.n}_m1={self.m1}_m2={self.m2}_p={self.p:.2f}"
