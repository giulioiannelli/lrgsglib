import networkx as nx

from ....config.const import *
from ..funcs import *
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX
from .generators_dgm import dorogovtsev_goltsev_mendes_graph_FastPatch


class DGMgraphNX(SignedGraphNX):
    """Dorogovtsev-Goltsev-Mendes graph as a :class:`SignedGraphNX`."""

    def __init__(
        self,
        n: int,
        *,
        with_positions: bool = False,
        sgpathn: str = DGM_SGPATH,
        stdFnameSFFX: str = DGM_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.n = n
        self.with_positions = with_positions
        self.only_const_mode = only_const_mode
        self.sgpathn = DGM_PHTABB if not sgpathn else sgpathn
        self.__init_stdFname__(stdFnameSFFX)
        if not only_const_mode:
            self.__init_network__()
        else:
            self.syshape = 3 * 2**n - 2
            self.syshapePth = f"N={self.syshape}"
            self.G = nx.Graph()
        super(DGMgraphNX, self).__init__(self.G, **kwargs)

    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.std_fname = DGM_PHTABB + SFFX

    def __init_network__(self) -> None:
        self.H = dorogovtsev_goltsev_mendes_graph_FastPatch(self.n)
        self.G = nx.convert_node_labels_to_integers(self.H)
        if self.with_positions:
            pos = nx.spring_layout(self.H)
            nx.set_node_attributes(self.H, pos, "pos")
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={self.syshape}"

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the DGM graph."""
        return int(self.syshape)

    def get_central_edge(self, on_g: str = DGM_ONREP):
        return (0, 1)
