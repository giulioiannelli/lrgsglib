from ..common import *
from ..funcs import *
from ..SignedGraph.SignedGraph import SignedGraph
from .generators import *
#
class Lattice2D(SignedGraph):
    eta_c = 1.128
    #
    def __init__(
        self,
        side1: int = L2D_SIDE1,
        geo: str = L2D_GEO,
        side2: int = L2D_SIDE2,
        pbc: bool = L2D_PBC,
        fbc_val: float = L2D_FBCV,
        stdFnameSFFX: str = L2D_STDFN,
        sgpathn: str = L2D_SGPATH,
        with_positions: bool = L2D_WITH_POS,
        bend_positions: bool = L2D_BEND_POS,
        prew: float = L2D_PREW,
        only_const_mode: bool = L2D_ONLY_CONST_MODE,
        **kwargs,
    ) -> None:
        self._verify_pflip(kwargs.get('pflip', 0.))
        self.only_const_mode = only_const_mode
        self.__init_side__(side1, side2)
        self.pbc = pbc
        self.prew = prew
        self.fbc_val = fbc_val
        self.__init_geo__(geo)
        #
        #
        self.__init_stdFname__(stdFnameSFFX)
        #
        self.sgpathn = pth_join(sgpathn, L2D_PATH_DICT[self.geo]) if sgpathn else L2D_PATH_DICT[self.geo]
        self.with_positions = with_positions
        self.bend_positions = bend_positions
        #
        if not only_const_mode:
            self.__init_lattice__()
        else:
            self.G = nx.Graph()
        super(Lattice2D, self).__init__(self.G, **kwargs)
    #
    def __init_side__(self, side1: int, side2: int) -> None:
        """
        Initialize the sides of the lattice, ensuring they are positive integers.

        Parameters
        ----------
        side1 : int
            The first side of the lattice.
        side2 : int
            The second side of the lattice.

        Raises
        ------
        ValueError
            If either side1 or side2 is not a positive integer.
        """
        if not is_positive_int(side1):
            raise ValueError("side1 must be a positive integer.")
        if side2 and not is_positive_int(side2):
            raise ValueError("side2 must be a positive integer.")
        if side2:
            self.side1, self.side2 = (side2, side1) if side2 > side1 else (side1, side2)
        else:
            self.side1 = side1
        #
    #
    def __init_geo__(self, geo: str) -> None:
        self.geo = geo
        if self.prew > 0.:
            self.geo = geo + '_sw'
        if geo not in L2D_GEO_LIST:
            if geo not in L2D_GEO_SHRT_LIST:
                warnings.warn(L2D_WARNMSG_GEO, Lattice2DWarning)
                self.geo = L2D_GEO
            else:
                self.geo = L2D_SHRT_GEO_DICT[self.geo]
        if not hasattr(self, 'side2'):
            if self.geo == 'hexagonal':
                self.side2 = self.side1
                self.side1 = adjust_to_even(self.side1/np.sqrt(3))
                if (self.side1 % 2 or self.side2 % 2) and self.pbc:
                    raise ValueError(L2D_ERRMSG_GEO)
            else:
                self.side2 = self.side1
            # if (self.side1 % 2 or self.side2 % 2) and self.pbc:
            #     raise ValueError(DEFLattice2D_geoerrmsg)
                #
        # For rhomb_octagonal graphs each rhomb contributes four nodes
        self.node_multiplier = 4 if self.geo.startswith(L2D_SHRT_GEO_DICT['oct_sqr']) else 1
        total_nodes = self.node_multiplier * self.side1 * self.side2
        
        if self.side1 == self.side2:
            self.syshapePth = f"N={total_nodes}"
        elif self.side1 > self.side2:
            self.syshapePth = f"L1={self.side1}_L2={self.side2}"
        elif self.side2 > self.side1:
            self.syshapePth = f"L1={self.side2}_L2={self.side1}"
        if self.prew > 0.:
            self.syshapePth = self.syshapePth + f"_prew={self.prew:.3g}"
        
        # Set syshape to reflect actual number of nodes
        if self.node_multiplier > 1:
            # For oct_sqr: each side gets multiplied by 2 (since 2*2 = 4 nodes per rhomb)
            self.syshape = (2 * self.side1, 2 * self.side2)
        else:
            self.syshape = (self.side1, self.side2)
        #
        self.p_c = L2D_P_C_DICT[self.geo]
        self.r_c = np.sqrt(self.eta_c/(np.pi*self.p_c))
    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.std_fname = L2D_GEO_SHRT_DICT[self.geo] + SFFX
    #
    def __init_lattice__(self) -> None:
        #
        if self.geo.startswith(L2D_SHRT_GEO_DICT['tri']):
            if self.prew == 0.:
                nxfunc = triangular_lattice_graph_FastPatch
            else:
                nxfunc = compose(triangular_lattice_graph_FastPatch, 
                                 rewire_edges_optimized, 
                                 g_kwargs={'prew': self.prew})
            self.z = 6
        elif self.geo.startswith(L2D_SHRT_GEO_DICT['sqr']):
            if self.prew == 0.:
                nxfunc = squared_lattice_graph_FastPatch
            else:
                nxfunc = compose(squared_lattice_graph_FastPatch, 
                                 rewire_edges_optimized, 
                                 g_kwargs={'prew': self.prew})
            self.z = 4
        elif self.geo == L2D_SHRT_GEO_DICT['hex']:
            self.z = 3
            nxfunc = hexagonal_lattice_graph_FastPatch
        elif self.geo.startswith(L2D_SHRT_GEO_DICT['oct_sqr']):
            self.z = 3
            nxfunc = rhomb_octagonal_graph_FastPatch

        #
        self.H = nxfunc(
            self.side1, 
            self.side2, 
            periodic=self.pbc, 
            with_positions=self.with_positions,
            bend_positions=self.bend_positions
        )
        self.G = nx.convert_node_labels_to_integers(self.H)

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the lattice."""
        return int(self.node_multiplier * self.side1 * self.side2)
    #
    def degree_check(self, degree):
        return np.where(np.array(list(map(lambda x: x[1], 
                                          list(self.G.degree())))) != degree)
    #
    def get_central_edge(self, on_g: str = L2D_ONREP):
        cnode = (self.side1//2-1, self.side2//2)
        cnode_t = (self.side1//2, self.side2//2)
        if self.geo == 'triangular':
            cnode = (self.side2//2, self.side1//2-1)
            cnode_t = (self.side2//2, self.side1//2)
        edge_t = (cnode, cnode_t)
        if not self.H.has_edge(*edge_t):
            if self.geo =='hexagonal':
                cnode = cnode_t
                cnode_t = (self.side1//2+1, self.side2//2)
                edge_t = (cnode, cnode_t)
        if on_g == 'G':
            return self.map_edge['G']['H'][edge_t]
        elif on_g == 'H':
            return edge_t
    #
    class nwContainer(dict):
        def __init__(self, l: SignedGraph, iterable=[], constant=None, 
                     **kwargs):
            super().__init__(**kwargs)
            self.update((key, constant) for key in iterable)
            self.l = l
            self.rd = self.l.graph_reprs
            self.rNodeFlip = {g: random.sample(list(self.l.gr[g].nodes()), 
                                    self.l.nflip
                                ) for g in self.rd}
            #
            self.centedge = {g: self.l.get_central_edge(g) 
                             for g in self.rd}
            self['single'] = {g: [self.centedge[g]] for g in self.rd}
            self['singleZERR'] = {g: self.get_links_ZERR(
                self.centedge[g][0], g, self.l.geo) for g in self.rd}
            self['singleXERR'] = {g: self.get_links_XERR(
                self.centedge[g][0], g) for g in self.rd}
            self['rand'] = {g: [e for e in self.l.fleset[g]] 
                            for g in self.rd}
            self['randZERR'] = {g: self.get_rand_pattern('ZERR', on_g=g) 
                                   for g in self.rd}
            self['randXERR'] = {g: self.get_rand_pattern('XERR', on_g=g) 
                                   for g in self.rd}
        #
        def get_links_XERR(self, node: Any, on_g: str = L2D_ONREP):
            return [(node, nn) for nn in self.l.get_graph_neighbors(node, on_g)]
        #
        def get_links_ZERR(self, node: Any, on_g: str = L2D_ONREP, 
                           geometry: str = L2D_GEO):
            dd = {'triangular': self.get_links_triangle,
             'squared': self.get_links_square,
             'hexagonal': self.get_links_hexagon}
            return dd[geometry](node, on_g)
        #
        def get_links_triangle(self, node: Any, on_g = L2D_ONREP):
            node2 = list(self.l.get_graph_neighbors(node, on_g))[0]
            common_neighbors = list(nx.common_neighbors(
                self.l.gr[on_g], node, node2))
            try:
                node3 = common_neighbors[0]
                links = [(node, node2), (node2, node3), (node, node3)]
            except IndexError:
                links = [(node, node2)]
            return links
        #
        def get_links_square(self, node: Any, on_g = L2D_ONREP):
            g = self.l.gr[on_g]
            neighbors = list(g.neighbors(node))
            for i in range(1, len(neighbors)):
                first_neighbor = neighbors[0]  # Always the first neighbor
                second_neighbor = neighbors[i]  # Iterating through the rest
                # Find common neighbors excluding the start_node
                common_neighbors = set(g.neighbors(first_neighbor)) \
                    & set(g.neighbors(second_neighbor))
                common_neighbors.discard(node)
                # If there is a common neighbor, we found a square
                if common_neighbors:
                    common_neighbor = common_neighbors.pop()
                    # Now, extract the links forming the square
                    links = [(node, first_neighbor),
                            (node, second_neighbor),
                            (first_neighbor, common_neighbor),
                            (second_neighbor, common_neighbor)]
                    return links
            links = [(node, first_neighbor),
                     (node, second_neighbor)]
            return links
        #
        def get_links_hexagon(self, node: int, 
                              on_g: str = L2D_ONREP):
            graph = self.l.gr[on_g]
            nodes_in_cycle = [node]
            node_nn = list(self.l.gr[on_g].neighbors(node))

            samp_node_nn_1 = node_nn[0]
            node_nn.remove(samp_node_nn_1)
            # #
            node_nn_1 = list(graph.neighbors(samp_node_nn_1))
            node_nn_1.remove(node)
            samp_node_nn_2 = node_nn_1[0]
            node_nn_1.remove(samp_node_nn_2)
            #
            flag = True
            for nn in node_nn:
                if flag:
                    common_nn = None
                    node_nn_1b = list(graph.neighbors(nn))
                    node_nn_1b.remove(node)
                    for i in node_nn_1b:
                        common_neighs = list(
                            nx.common_neighbors(graph, samp_node_nn_2, i))
                        if common_neighs != []:
                            nodes_in_cycle.extend([nn, i, common_neighs[0]])
                            flag = False
                            break
            try:
                nodes_in_cycle.extend([samp_node_nn_1, samp_node_nn_2])
            except IndexError:
                pass
            subH = graph.subgraph(nodes_in_cycle)
            links = [tuple(sorted(edge)) for edge in subH.edges()]
            return links
        #
        def get_rand_pattern(self, mode: str, on_g: str = L2D_ONREP):
            match mode:
                case "XERR":
                    if COUNT_XERR_PATTERNS:
                        patternList = [k for i in self.rNodeFlip[on_g] 
                                    for k in self.get_links_XERR(i, on_g)]
                    else:
                        tmplst = self.rNodeFlip[on_g]
                        grph = self.l.gr[on_g]
                        _ = 0
                        patternList = []
                        while _ < len(tmplst):
                            leval = [all([nnn['weight'] == -1 
                                        for nnn in grph[nn].values()])
                                        for nn in grph.neighbors(tmplst[_])]
                            if any(leval):
                                tmplst.pop(_)  # Removing the element
                            else:
                                glXERR = self.get_links_XERR(tmplst[_], 
                                                             on_g)
                                patternList.extend([k for k in glXERR])
                                _ += 1
                case "ZERR":
                    match self.l.geo:
                        case 'squared':
                            patternList = [k for i in self.rNodeFlip[on_g] 
                                for k in self.get_links_square(i, on_g)]
                        case "hexagonal":
                            patternList = [k for i in self.rNodeFlip[on_g]
                                for k in self.get_links_hexagon(i, on_g)]    
                        case "triangular":
                            patternList = [k for i in self.rNodeFlip[on_g] 
                               for k in self.get_links_triangle(i, on_g)]
            return list(set(patternList))
        #
        def get_links_rball(self, R: int = 1, center: Any = None, 
                            on_g: str = L2D_ONREP):
            graph = self.l.gr[on_g]
            if not center:
                center = self.centedge[on_g][0]
            neighs_to_flip = get_neighbors_at_distance(graph, center, R)
            links = {(node, neighbor) for node in neighs_to_flip 
                     for neighbor in graph.neighbors(node)}
            return links
    #     # #
    # def make_animation(self, fig, ax, frames):

    #     cv0 = frames[0].reshape(self.syshape)
    #     im = ax.imshow(cv0)  # Here make an AxesImage rather than contour
    #     _, _, cbar = imshow_colorbar_caxdivider(im, ax)
    #     # tx = ax.set_title('Frame 0')

    #     def animate(i):
    #         arr = frames[i].reshape(self.syshape)
    #         vmax = np.max(arr)
    #         vmin = np.min(arr)
    #         im.set_data(arr)
    #         cbar.mappable.set_clim(vmin, vmax)
    #         # tx.set_text('Frame {0}'.format(i))
    #         # In this version you don't have to do anything to the colorbar,
    #         # it updates itself when the mappable it watches (im) changes

    #     return animate