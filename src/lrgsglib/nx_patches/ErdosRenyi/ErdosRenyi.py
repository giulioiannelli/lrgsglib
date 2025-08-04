from ..common import *
from ..funcs import *
from ..SignedGraph.SignedGraph import SignedGraph


class ErdosRenyi(SignedGraph):
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
        super(ErdosRenyi, self).__init__(self.G, **kwargs)

    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.std_fname = ER_PHTABB + SFFX

    #
    def __init_network__(self, n, p):
        G = nx.erdos_renyi_graph(n, p)
        CC = nx.connected_components(G)
        GC = max(CC, key=len)
        self.G = G.subgraph(GC).copy()
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={n}_p={p:.3g}"

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the Erdos-Rényi graph."""
        return int(self.syshape)

    class nwContainer(dict):
        def __init__(
            self, er: SignedGraph, iterable=[], constant=None, **kwargs
        ):
            super().__init__(**kwargs)
            self.update((key, constant) for key in iterable)
            self.er = er
            self.rd = self.er.graph_reprs
            self.rNodeFlip = {
                g: random.sample(list(self.er.nodes_in(g)), self.er.nflip)
                for g in self.rd
            }
            self["rand"] = {g: [e for e in self.er.fleset[g]] for g in self.rd}
            self["randXERR"] = {
                g: self.get_rand_pattern("XERR", on_g=g) for g in self.rd
            }

        #
        def get_links_XERR(self, node: Any, on_g: str = ER_ONREP):
            return [
                (node, nn) for nn in self.er.get_graph_neighbors(node, on_g)
            ]

        #
        def get_links_ZERR(self, node: Any, on_g: str = ER_ONREP):
            smallest_cycle = get_smallest_cycle_graph_node(self.rd[on_g], node)
            edges = []
            if smallest_cycle:
                edges = [
                    (smallest_cycle[i], smallest_cycle[i + 1])
                    for i in range(len(smallest_cycle) - 1)
                ]
            return edges

        #
        def get_rand_pattern(self, mode: str, on_g: str = ER_ONREP):
            match mode:
                case "XERR":
                    if COUNT_XERR_PATTERNS:
                        patternList = [
                            k
                            for i in self.rNodeFlip[on_g]
                            for k in self.get_links_XERR(i, on_g)
                        ]
                    else:
                        tmplst = self.rNodeFlip[on_g]
                        grph = self.er.gr[on_g]
                        _ = 0
                        patternList = []
                        while _ < len(tmplst):
                            leval = [
                                all(
                                    [
                                        nnn["weight"] == -1
                                        for nnn in grph[nn].values()
                                    ]
                                )
                                for nn in grph.neighbors(tmplst[_])
                            ]
                            if any(leval):
                                tmplst.pop(_)  # Removing the element
                            else:
                                glXERR = self.get_links_XERR(tmplst[_], on_g)
                                patternList.extend([k for k in glXERR])
                                _ += 1
                case "ZERR":
                    tmplst = self.rNodeFlip[on_g]
                    grph = self.er.gr[on_g]
                    _ = 0
                    patternList = [
                        links
                        for node in tmplst
                        if (links := self.get_links_ZERR(node, on_g))
                    ]
                    # for node in tmplst:
                    #     patternList.extend(self.get_links_XERR(node, on_g))
            return patternList
