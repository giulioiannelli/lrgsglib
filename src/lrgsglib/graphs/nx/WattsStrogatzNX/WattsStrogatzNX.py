import networkx as nx

from ....config.const import *
from ..funcs import *
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


class WattsStrogatzNX(SignedGraphNX):
    """
    Signed Watts-Strogatz small-world graph generator.

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Each node is connected to k nearest neighbors in ring topology.
    p : float
        Probability of rewiring each edge.
    sgpathn : str, default WS_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default WS_STDFN
        Filename suffix used in on-disk exports.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraphNX`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Notes
    -----
    The graph is generated using NetworkX's ``connected_watts_strogatz_graph``
    which guarantees a connected graph.

    ``pflip`` defines how many edges are marked for sign flips, but weights
    are not set negative until you call ``flip_random_fract_edges`` or
    ``flip_sel_edges``, unless weights already exist on the input graph.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import WattsStrogatzNX
    >>> ws = WattsStrogatzNX(n=100, k=4, p=0.3, pflip=0.2, seed=1)
    >>> ws.flip_random_fract_edges()  # apply sign flips
    """
    def __init__(
        self,
        n: int,
        k: int,
        p: float,
        sgpathn: str = WS_SGPATH,
        stdFnameSFFX: str = WS_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ):

        self.sgpathn = WS_PHTABB if not sgpathn else sgpathn
        self.only_const_mode = only_const_mode
        self.__init_stdFname__(stdFnameSFFX)
        if not only_const_mode:
            self.__init_network__(n, k, p)
        else:
            self.n = n
            self.k = k
            self.p = p
            self.syshape = n
            self.syshapePth = f"N={n}_p={p:.3g}"
            self.G = nx.Graph()
        super(WattsStrogatzNX, self).__init__(self.G, **kwargs)

    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.std_fname = WS_PHTABB + SFFX

    #
    def __init_network__(self, n, k, p):
        self.G = nx.connected_watts_strogatz_graph(n, k, p)
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={n}_p={p:.3g}"

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the Watts-Strogatz graph."""
        return int(self.syshape)
