"""
RandomGeometric: Random geometric graph with spatial constraints.

Nodes are placed uniformly in a unit hypercube and connected if their
distance is less than a threshold radius. This creates graphs with
spatial clustering and percolation transitions.

Example
-------
>>> from lrgsglib.nx_patches.random import RandomGeometric
>>> rgg = RandomGeometric(n=200, radius=0.15, dim=2, pflip=0.1, seed=42)
>>> rgg.flip_random_fract_edges()
"""

import networkx as nx

from .random_graph import RandomGraph


# Constants
RGG_PHTABB = "rgg"
RGG_SGPATH = ""
RGG_STDFN = ""


class RandomGeometric(RandomGraph):
    """
    Signed random geometric graph.

    Nodes are placed uniformly at random in a unit hypercube [0,1]^dim.
    Two nodes are connected if their Euclidean distance is less than
    the given radius.

    Parameters
    ----------
    n : int
        Number of nodes.
    radius : float
        Connection threshold distance.
    dim : int, default 2
        Dimension of the embedding space.
    sgpathn : str, default RGG_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default RGG_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraph (pflip, seed, etc.).

    Attributes
    ----------
    radius : float
        Connection threshold.
    dim : int
        Embedding dimension.

    Notes
    -----
    Random geometric graphs are useful for:
    - Modeling wireless networks and sensor networks
    - Studying percolation transitions
    - Spatial network analysis
    - Graphs with natural clustering

    The critical radius for connectivity in 2D is approximately
    r_c ~ sqrt(log(n) / (pi * n)).

    Node positions are stored as 'pos' attribute if with_positions=True.

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import RandomGeometric
    >>> import numpy as np
    >>> # Create RGG near percolation threshold
    >>> n = 500
    >>> r_c = np.sqrt(np.log(n) / (np.pi * n))
    >>> rgg = RandomGeometric(n=n, radius=1.5*r_c, dim=2, seed=42)
    >>> print(f"Nodes: {rgg.N}, Connected: {rgg.N == n}")
    """

    def __init__(
        self,
        n: int,
        radius: float,
        dim: int = 2,
        sgpathn: str = RGG_SGPATH,
        stdFnameSFFX: str = RGG_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius}")
        if dim < 1:
            raise ValueError(f"dim must be >= 1, got {dim}")

        self.radius = radius
        self.dim = dim
        self._init_std_fname(stdFnameSFFX)
        sgpath = RGG_PHTABB if not sgpathn else sgpathn

        # RGG stores positions automatically
        super().__init__(
            n=n,
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=True,  # May be disconnected below percolation
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        self.std_fname = RGG_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate random geometric graph."""
        return nx.random_geometric_graph(self.n, self.radius, dim=self.dim)

    def _compute_syshapePth(self) -> str:
        return f"N={self.n}_r={self.radius:.3g}_d={self.dim}"

    def get_critical_radius(self) -> float:
        """
        Return the approximate critical radius for connectivity in 2D.

        Returns
        -------
        float
            Critical radius r_c ~ sqrt(log(n) / (pi * n)).
        """
        import numpy as np
        if self.dim != 2:
            raise NotImplementedError("Critical radius formula only for dim=2")
        return float(np.sqrt(np.log(self.n) / (np.pi * self.n)))

    def get_positions(self) -> dict:
        """
        Return node positions from the graph.

        Returns
        -------
        dict
            Dictionary mapping node to (x, y, ...) position tuple.
        """
        return nx.get_node_attributes(self.G, 'pos')
