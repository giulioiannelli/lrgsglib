"""
HolmeKimNX: Holme-Kim model: Barabasi-Albert with clustering.

Extends BA by adding a "triad formation" step. After connecting to a
node via preferential attachment, with probability p the new node also
connects to a neighbor of that node, forming triangles.

Example
-------
>>> from lrgsglib.graphs.nx import HolmeKimNX
>>> g = HolmeKimNX(n=500, m=3, p=0.8, seed=42)
>>> import networkx as nx
>>> print(f"Clustering: {nx.average_clustering(g.G):.3f}")
"""

import networkx as nx

from ..RandomGraphNX.RandomGraphNX import RandomGraphNX


# Constants
PA_PHTABB = "pa"
PA_SGPATH = ""
PA_STDFN = ""


class HolmeKimNX(RandomGraphNX):
    """
    Holme-Kim model: Barabasi-Albert with clustering.

    Extends BA by adding a "triad formation" step. After connecting to a
    node via preferential attachment, with probability p the new node also
    connects to a neighbor of that node, forming triangles.

    Parameters
    ----------
    n : int
        Number of nodes.
    m : int
        Number of edges to attach from each new node.
    p : float, default 0.5
        Probability of triad formation step.
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
    - p=0: Standard BA (no extra clustering)
    - p=1: Maximum clustering (always try to form triangles)

    This model produces scale-free networks with higher clustering
    coefficient than standard BA, more realistic for social networks.

    References
    ----------
    Holme, P., & Kim, B. J. (2002). Growing scale-free networks with
    tunable clustering. Physical Review E, 65(2), 026107.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import HolmeKimNX
    >>> # High clustering scale-free
    >>> g = HolmeKimNX(n=500, m=3, p=0.8, seed=42)
    >>> import networkx as nx
    >>> print(f"Clustering: {nx.average_clustering(g.G):.3f}")
    """

    def __init__(
        self,
        n: int,
        m: int,
        p: float = 0.5,
        sgpathn: str = PA_SGPATH,
        stdFnameSFFX: str = PA_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if m < 1:
            raise ValueError(f"m must be >= 1, got {m}")
        if not 0 <= p <= 1:
            raise ValueError(f"p must be in [0, 1], got {p}")
        if m >= n:
            raise ValueError(f"m must be < n, got m={m}, n={n}")

        self.m = m
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
        self.std_fname = PA_PHTABB + "_hk" + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate Holme-Kim graph."""
        seed = getattr(self, '_rng_seed', None)
        return nx.powerlaw_cluster_graph(self.n, self.m, self.p, seed=seed)

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"N={self.n}_m={self.m}_p={self.p:.2f}"

    def get_clustering_coefficient(self) -> float:
        """
        Return the average clustering coefficient.

        Returns
        -------
        float
            Average clustering coefficient.
        """
        return nx.average_clustering(self.G)
