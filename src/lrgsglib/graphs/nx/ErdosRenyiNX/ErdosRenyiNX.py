import networkx as nx

from ....config.const import *
from ..funcs import *
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


class ErdosRenyiNX(SignedGraphNX):
    """
    Signed Erdos-Renyi graph generator (largest connected component).

    Parameters
    ----------
    n : int
        Number of nodes in the initial Erdos-Renyi graph.
    p : float
        Edge probability for the Erdos-Renyi model.
    sgpathn : str, default ER_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default ER_STDFN
        Filename suffix used in on-disk exports.
    only_const_mode : bool, default False
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraphNX`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Notes
    -----
    The graph is generated with NetworkX and then reduced to its largest
    connected component before initializing ``SignedGraphNX``. This ensures
    that most downstream algorithms operate on a connected graph.

    ``pflip`` defines how many edges are marked for sign flips, but weights
    are not set negative until you call ``flip_random_fract_edges`` or
    ``flip_sel_edges``, unless weights already exist on the input graph.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import ErdosRenyiNX
    >>> er = ErdosRenyiNX(n=200, p=0.05, pflip=0.2, seed=1)
    >>> er.flip_random_fract_edges()  # apply sign flips
    """
    def __init__(
        self,
        n: int,
        p: float,
        sgpathn: str = ER_SGPATH,
        stdFnameSFFX: str = ER_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ):
        self.only_const_mode = only_const_mode
        self.sgpathn = ER_PHTABB if not sgpathn else sgpathn
        self.__init_stdFname__(stdFnameSFFX)
        self.n = n
        self.p = p
        if not only_const_mode:
            self.__init_network__(n, p)
        else:
            self.syshape = n
            self.syshapePth = f"N={n}_p={p:.3g}"
            self.G = nx.Graph()
        super(ErdosRenyiNX, self).__init__(self.G, **kwargs)

    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.std_fname = ER_PHTABB + SFFX

    #
    def __init_network__(self, n, p):
        G = nx.erdos_renyi_graph(n, p)
        CC = nx.connected_components(G)
        GC = max(CC, key=len)
        # Relabel the giant component to contiguous 0..N-1 integers. Dropping
        # isolated nodes otherwise leaves gaps in the node labels, which breaks
        # every positional consumer (e.g. the statsys sweepers index the state
        # array by node id); contiguous labels are the invariant all other
        # graph types already satisfy.
        self.G = nx.convert_node_labels_to_integers(G.subgraph(GC).copy())
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={n}_p={p:.3g}"

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the Erdos-Renyi graph."""
        return int(self.syshape)

    # nwDict patterns: inherit the engine-neutral ``NwContainer`` from
    # ``SignedGraphNX``. The old bespoke container only built ``rand`` /
    # ``randXERR`` (its ``get_links_ZERR`` was dead code); the shared container
    # builds those identically *and* supplies the structured ``single*`` /
    # ``*ZERR`` patterns on demand (via the base hub ``get_central_edge``), so
    # disorder supports like ``Disorder('singleXERR')`` now work on ER too.
