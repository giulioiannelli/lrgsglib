from ..common import *
from ..funcs import *
from ..SignedGraph.SignedGraph import SignedGraph
from .generators import dorogovtsev_goltsev_mendes_graph_FastPatch


class DGMgraph(SignedGraph):
    """Dorogovtsev-Goltsev-Mendes graph as a :class:`SignedGraph`."""

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
        super(DGMgraph, self).__init__(self.G, **kwargs)

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

    class nwContainer(dict):
        def __init__(self, dgm: SignedGraph, iterable=[], constant=None, **kwargs):
            super().__init__(**kwargs)
            self.update((key, constant) for key in iterable)
            self.dgm = dgm
            self.rd = self.dgm.graph_reprs
            self.rNodeFlip = {
                g: random.sample(list(self.dgm.nodes_in(g)), self.dgm.nflip)
                for g in self.rd
            }
            self['rand'] = {g: [e for e in self.dgm.fleset[g]] for g in self.rd}
            self['randXERR'] = {g: self.get_rand_pattern('XERR', on_g=g) for g in self.rd}

        def get_links_XERR(self, node: Any, on_g: str = DGM_ONREP):
            return [(node, nn) for nn in self.dgm.get_graph_neighbors(node, on_g)]

        def get_rand_pattern(self, mode: str, on_g: str = DGM_ONREP):
            match mode:
                case 'XERR':
                    if COUNT_XERR_PATTERNS:
                        patternList = [k for i in self.rNodeFlip[on_g]
                                        for k in self.get_links_XERR(i, on_g)]
                    else:
                        tmplst = self.rNodeFlip[on_g]
                        grph = self.dgm.gr[on_g]
                        _ = 0
                        patternList = []
                        while _ < len(tmplst):
                            leval = [all([nnn['weight'] == -1 for nnn in grph[nn].values()])
                                     for nn in grph.neighbors(tmplst[_])]
                            if any(leval):
                                tmplst.pop(_)
                            else:
                                glXERR = self.get_links_XERR(tmplst[_], on_g)
                                patternList.extend([k for k in glXERR])
                                _ += 1
            return list(set(patternList))
