"""
BipartiteFromDegreeSequenceNX: Signed bipartite graph with prescribed degree sequences.

Creates a bipartite graph where top and bottom nodes have specified degree sequences.

Example
-------
>>> from lrgsglib.graphs.nx import BipartiteFromDegreeSequenceNX
>>> # Create bipartite with specified degrees
>>> top_deg = [3, 3, 2, 2]  # sum = 10
>>> bot_deg = [2, 2, 2, 2, 2]  # sum = 10
>>> bg = BipartiteFromDegreeSequenceNX(top_deg, bot_deg, seed=42)
"""

from typing import Sequence

import networkx as nx
from networkx.algorithms import bipartite

from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


# Constants
BP_PHTABB = "bp"
BP_SGPATH = ""
BP_STDFN = ""


class BipartiteFromDegreeSequenceNX(SignedGraphNX):
    """
    Signed bipartite graph with prescribed degree sequences.

    Creates a bipartite graph where top and bottom nodes have
    specified degree sequences.

    Parameters
    ----------
    top_degrees : Sequence[int]
        Degree sequence for top nodes.
    bottom_degrees : Sequence[int]
        Degree sequence for bottom nodes.
    sgpathn : str, default BP_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default BP_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraphNX (pflip, seed, etc.).

    Notes
    -----
    The sum of top_degrees must equal the sum of bottom_degrees
    (total number of edges from each side must match).

    Examples
    --------
    >>> from lrgsglib.graphs.nx import BipartiteFromDegreeSequenceNX
    >>> # Create bipartite with specified degrees
    >>> top_deg = [3, 3, 2, 2]  # sum = 10
    >>> bot_deg = [2, 2, 2, 2, 2]  # sum = 10
    >>> bg = BipartiteFromDegreeSequenceNX(top_deg, bot_deg, seed=42)
    """

    def __init__(
        self,
        top_degrees: Sequence[int],
        bottom_degrees: Sequence[int],
        sgpathn: str = BP_SGPATH,
        stdFnameSFFX: str = BP_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.top_degrees = list(top_degrees)
        self.bottom_degrees = list(bottom_degrees)

        if sum(self.top_degrees) != sum(self.bottom_degrees):
            raise ValueError(
                f"Sum of degree sequences must be equal: "
                f"sum(top)={sum(self.top_degrees)} != sum(bottom)={sum(self.bottom_degrees)}"
            )

        self.n1 = len(top_degrees)
        self.n2 = len(bottom_degrees)
        self.only_const_mode = only_const_mode
        self._init_std_fname(stdFnameSFFX)
        self.sgpathn = BP_PHTABB if not sgpathn else sgpathn
        self.syshape = (self.n1, self.n2)
        self.syshapePth = self._compute_syshapePth()

        if not only_const_mode:
            self._init_bipartite()
        else:
            self.G = nx.Graph()
            self.top_nodes = set()
            self.bottom_nodes = set()

        super().__init__(self.G, **kwargs)

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = BP_PHTABB + "_deg" + suffix

    def _init_bipartite(self) -> None:
        """Build the bipartite graph from degree sequences."""
        # Use Havel-Hakimi-style algorithm
        self.G = bipartite.configuration_model(
            self.top_degrees, self.bottom_degrees, create_using=nx.Graph()
        )
        # Remove parallel edges and self-loops
        self.G = nx.Graph(self.G)

        # Store node sets
        self.top_nodes = {n for n, d in self.G.nodes(data=True) if d["bipartite"] == 0}
        self.bottom_nodes = {
            n for n, d in self.G.nodes(data=True) if d["bipartite"] == 1
        }

    def _compute_syshapePth(self) -> str:
        """Compute system shape path string."""
        return f"n1={self.n1}_n2={self.n2}_deg"

    def get_expected_num_nodes(self) -> int:
        """Return total number of nodes."""
        return self.n1 + self.n2

    def get_expected_num_edges(self) -> int:
        """Return expected number of edges."""
        return sum(self.top_degrees)
