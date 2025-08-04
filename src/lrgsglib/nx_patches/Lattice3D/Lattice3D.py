# from ..common import *
import random
import warnings
#
from typing import Union, Tuple, Any
from networkx import convert_node_labels_to_integers, set_node_attributes, Graph
#
from os.path import join as pth_join
#
from ...config.const import *
from ...config.errwar import Lattice2DWarning
import numpy as np
from ...utils.basic.geometry import project_3d_to_2d
from ...utils.basic.numeric import is_positive_int
from ...utils.basic.functions import compose
from ..funcs import LatticeND_graph_FastPatch, remove_edges
from ..SignedGraph.SignedGraph import SignedGraph
from .generators import *

class Lattice3D(SignedGraph):
    def __init__(
        self,
        dim: Union[int, Tuple[int, int, int]] = L3D_DIM,
        geo: str = L3D_GEO,
        pbc: bool = L3D_PBC,
        fbc_val: float = L3D_FBCV,
        stdFnameSFFX: str = L3D_STDFN,
        sgpathn: str = L3D_SGPATH,
        with_positions: bool = L3D_WITH_POS,
        pdil: float = L3D_PDIL,
        theta: float = L3D_THETA,
        phi: float = L3D_PHI,
        only_const_mode: bool = L3D_ONLY_CONST_MODE,
        **kwargs,
    ) -> None:
        self.only_const_mode = only_const_mode
        self.__init_dim__(dim)
        self.pbc = pbc
        self.pdil = pdil
        self.fbc_val = fbc_val
        self.geo = L3D_GEO_DICT[geo]
        self.theta = theta
        self.phi = phi
        self.__init_stdFname__(stdFnameSFFX)
        _ = L3D_PATH_DICT[self.geo]
        self.sgpathn = pth_join(sgpathn, L3D_PATH_DICT[self.geo]) if sgpathn else L3D_PATH_DICT[self.geo]
        self.with_positions = with_positions
        if not only_const_mode:
            self.__init_lattice__()
        else:
            if self.geo == L3D_GEO_BCC:
                self.node_multiplier = 2
            elif self.geo == L3D_GEO_FCC:
                self.node_multiplier = 4
            else:
                self.node_multiplier = 1
            self.syshape = self.dim
            total_nodes = self.node_multiplier * int(np.prod(self.dim))
            if all(x == self.dim[0] for x in self.dim):
                self.syshapePth = f"N={total_nodes}"
            else:
                dim_part = '_'.join([f"L{i}={side}" for i, side in enumerate(self.dim)])
                self.syshapePth = f"{dim_part}_N={total_nodes}"
            self.G = Graph()
        super(Lattice3D, self).__init__(self.G, **kwargs)
    #
    def __init_dim__(self, dim: Union[int, Tuple[int, int, int]]) -> None:
        if is_positive_int(dim):
            self.dim = (dim, dim, dim)
        elif (
            isinstance(dim, tuple)
            and len(dim) == 3
            and all(is_positive_int(x) for x in dim)
        ):
            self.dim = tuple(sorted(dim, reverse=True))
        else:
            raise ValueError("dim must be a positive integer or a tuple of 3 positive integers")

        self.dimL = list(self.dim)
    #
    def __init_geo__(self, geo: str) -> None:
        self.geo = geo
        if self.pdil > 0.:
            self.geo = geo + '_dil'
        if geo not in L3D_GEO_LIST:
            if geo not in L3D_GEO_SHRT_LIST:
                warnings.warn(L3D_WARNMSG_GEO, Lattice2DWarning)
                self.geo = L3D_GEO
            else:
                self.geo = L3D_SHRT_GEO_DICT[self.geo]
    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.std_fname = L3D_STDFN + SFFX
    #
    def __init_lattice__(self) -> None:
        if self.geo == L3D_GEO_SC:
            self.node_multiplier = 1
            if self.pdil == 0.:
                nxfunc = LatticeND_graph_FastPatch
            else:
                nxfunc = compose(
                    LatticeND_graph_FastPatch,
                    remove_edges,
                    g_kwargs={'pdil': self.pdil},
                )
        elif self.geo == L3D_GEO_BCC:
            self.node_multiplier = 2
            nxfunc = generate_bcc_lattice
        elif self.geo == L3D_GEO_FCC:
            self.node_multiplier = 4
            nxfunc = generate_fcc_lattice
        else:
            raise ValueError(f"Unsupported geometry '{self.geo}'.")
        self.syshape = self.dim
        total_nodes = self.node_multiplier * int(np.prod(self.dim))
        if all(x == self.dim[0] for x in self.dim):
            self.syshapePth = f"N={total_nodes}"
        else:
            dim_part = '_'.join([f"L{i}={side}" for i, side in enumerate(self.dim)])
            self.syshapePth = f"{dim_part}_N={total_nodes}"
        
        self.H = nxfunc(self.dim, periodic=self.pbc)
        self.G = convert_node_labels_to_integers(self.H)

        if self.with_positions:
            self._set_positions()

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the 3D lattice."""
        return int(self.node_multiplier * np.prod(self.dim))
    #
    def _set_positions(self):
        pos = {node: project_3d_to_2d(*node, self.theta, self.phi)
               for node in self.H.nodes()}
        set_node_attributes(self.H, pos, 'pos')


    def get_central_edge(self, on_g: str = L3D_ONREP):
        cnode = (self.dimL[0]//2-1, self.dimL[1]//2, self.dimL[2]//2)
        cnode_t = (self.dimL[0]//2, self.dimL[1]//2, self.dimL[2]//2)
        edge_t = (cnode, cnode_t)
        if on_g == 'H':
            return edge_t
        elif on_g == 'G':
            return self.map_edge['G']['H'][edge_t]
    class nwContainer(dict):
        def __init__(self, l: SignedGraph, iterable=[], constant=None, 
                    **kwargs):
            super().__init__(**kwargs)
            self.update((key, constant) for key in iterable)
            self.l = l
            self.rd = self.l.graph_reprs
            self.rNodeFlip = {g: random.sample(
                                    list(self.l.gr[g].nodes()), 
                                    self.l.nflip
                                ) for g in self.rd}
            #
            self.centedge = {g: self.l.get_central_edge(g) 
                            for g in self.rd}
            self['single'] = {g: [self.centedge[g]] for g in self.rd}
            self['singleXERR'] = {g: self.get_links_XERR(
                self.centedge[g][0], g) for g in self.rd}
            self['rand'] = {g: [e for e in self.l.fleset[g]] 
                        for g in self.rd}
            self['randXERR'] = {g: self.get_rand_pattern('XERR', on_g=g) 
                        for g in self.rd}
        #
        def get_links_XERR(self, node: Any, on_g: str = L3D_ONREP):
            return [(node, nn) for nn in self.l.get_graph_neighbors(node, on_g)]
        #
        def get_rand_pattern(self, mode: str, on_g: str = L3D_ONREP):
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
            return list(set(patternList))
        

