from ..common import *
from ..SignedGraph.SignedGraph import SignedGraph


class WattStrogatz(SignedGraph):
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
        super(WattStrogatz, self).__init__(self.G, **kwargs)

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
