import copy
import os
import logging
import random
import struct
import time
import warnings

import cupy as cp
import networkx as nx
import numpy as np
import pickle as pk

from networkx import Graph
from typing import Any, Union, List, Dict
from matplotlib.pyplot import get_cmap
from numpy.typing import NDArray
from scipy.sparse import spdiags
from scipy.sparse import identity as scsp_identity
from scipy.sparse.linalg import eigsh as scsp_eigsh

from ...config.const import *
from ...config.errwar import NflipError, NoClustError, SignedGraphWarning
from ...utils.basic import is_in_range, join_non_empty,\
    flip_to_positive_majority_adapted, bin_sign, flip_to_positive_majority,\
    generate_random_id, normalize_array
from ...utils.lrg import compute_ising_pairwise_energy
from ...utils.tools import NestedDict, ConditionalPartitioning
#
logger = logging.getLogger(__name__)
class SignedGraph:
    sgpathn = "signed_graph"
    #
    def __init__(self,
        G: Graph, 
        pflip: float = SG_PFLIP,
        on_g: str = SG_GRAPH_REPR,
        seed: int = None,
        path_data: Path = None,
        path_plot: Path = None,
        init_nw_dict: bool = SG_INIT_NW_DICT,
        init_weights_val: Union[float, dict] = SG_INIT_WVAL,
        export_mode: str = SG_EXPORT_M,
        make_dir_tree: bool = True,
        imported: bool = SG_IMPORT_ON,
        import_fname: str = '',
        import_mode: str = SG_LOAD_M,
    ):
        self.graph_representation_dictionary = {}
        self.map_node = {}
        self.map_edge = {}
        self.eset = {}
        self.fleset = {}
        self.lfeset = {}
        #
        self._verify_pflip(pflip)
        self.__init_randomness__(seed)
        #
        self.__init_paths__(path_data, path_plot, make_dir_tree)
        #
        self.on_g = on_g
        self.init_weights_val = init_weights_val
        self.load_g = imported
        self.init_nw_dict = init_nw_dict
        self.std_fname = join_non_empty('_', self.std_fname, self.peq_str)
        #
        self.G = self.__load_graph__(import_fname, import_mode) \
            if self.load_g else G
        self.__init_reprdict__()
        if not self.only_const_mode:
            self.__init_sgraph__()
            if self.init_nw_dict:
                if hasattr(self, 'nwContainer'):
                    self.nwDict = self.nwContainer(self)
                else:
                    raise AttributeError(SG_ERRMSG_NW_DICT)
            self.__make_graph_clustering_utility__()
    #
    @property
    def adj(self):
        return self.adjacency_matrix
    @property
    def degm(self):
        return self.degree_matrix
    @property
    def gr(self):
        return self.graph_representation_dictionary
    @property
    def gcl(self):
        return self.graph_clustering_utility
    @property
    def lap(self):
        return self.laplacian_matrix
    @property
    def sdeg(self):
        return self.signed_degree_matrix
    @property
    def slp(self):
        return self.signed_laplacian_matrix
    @property
    def N(self):
        return self.gr[self.on_g].number_of_nodes()
    @property
    def Ne(self):
        return self.gr[self.on_g].number_of_edges()
    @property
    def Ne_n(self):
        return len(self.fleset[self.on_g])
    @property
    def nodes_in(self, on_g: str = SG_GRAPH_REPR):
        return list(self.gr[on_g].nodes())

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes without building the graph."""
        if hasattr(self, "syshape") and isinstance(self.syshape, int):
            return self.syshape
        raise NotImplementedError(
            "Subclasses must implement `get_expected_num_nodes`." )
    #
    def __init_randomness__(self, seed: int = None) -> None:
        """
        Initialize the random seed of both the random and numpy libraries.
        If no seed is provided, it will use the current time in milliseconds
        plus the process ID to ensure randomness across different runs.

        Parameters
        ----------
        seed : int, optional
            The seed value to initialize the random number generators. If not
            provided, a seed is generated using the current time and process ID.
        Returns
        -------
        None
        """
        if LB_PFLIP < self.pflip and self.pflip < UB_PFLIP:
            self.seed = seed or ((int(time.time() * 1000) + os.getpid()) % (2**32 - 1))
        else:
            self.seed = 0
        random.seed(self.seed)
        try:
            cp.random.seed(self.seed)
        except Exception as e:
            warnings.warn(f"Could not set cupy random seed: {e}. "
                          "Cupy will not be used for random operations.")
            pass
        np.random.seed(self.seed)
        #
        self.rand_str = generate_random_id()
        #
        logger.info(f"Random seed set to {self.seed}")
    #
    def __init_paths__(
            self,
            path_data: Path = None, 
            path_plot: Path = None,
            make_dir_tree: bool = True,
            exist_ok: bool = True
    ) -> None:
        """
        Initialize paths for data, plots, and various subdirectories used by the
        SignedGraph class. If no paths are provided, default paths are used.

        Parameters
        ----------
        path_data : Path, optional
            The base path for data storage. Defaults to the global PATHDATA.
        path_plot : Path, optional
            The base path for plot storage. Defaults to the global PATHPLOT.
        make_dir_tree : bool, optional
            If True, creates the necessary directory structure.
        exist_ok : bool, optional
            If True, existing directories will not raise an error.

        Returns
        -------
        dict
            A dictionary containing the initialized paths.
        """
        #
        self.path_data = path_data or PATHDATA
        self.path_plot = path_plot or PATHPLOT
        self.path_sgdata = self.path_data / Path(self.sgpathn)
        #
        self.subpath_list = []
        for p in PATHN_LIST:
            pfname = Path(p, self.syshapePth)
            setattr(self, f"path_{p}", self.path_sgdata / pfname)
            self.subpath_list.append(getattr(self, f"path_{p}"))
        #
        if make_dir_tree: self.__make_dirs__(exist_ok)
    #
    def __make_dirs__(self, exist_ok: bool = True):
        """
        Create the necessary directories for the SignedGraph class. The 
        directories are specified in the `subpath_list` attribute.

        Parameters
        ----------
        exist_ok : bool, optional
            If True, existing directories will not raise an error. 
        """
        for _ in self.subpath_list: 
            os.makedirs(_, exist_ok=exist_ok)
    #
    def __make_graph_clustering_utility__(self):
        self.graph_clustering_utility = NestedDict()
    #
    def __init_reprdict__(self):
        self.gr[self.on_g] = getattr(self, self.on_g)
        for grd in SG_LIST_REPR:
            try:
                self.gr[grd] = getattr(self, grd)
            except AttributeError:
                pass
        self.graph_reprs = list(self.gr.keys())

    def __init_loaded_graph__(
        self,
        path_data: Path = None, 
        path_plot: Path = None,
        on_g: str = SG_GRAPH_REPR
    ) -> None:
        self.graph_representation_dictionary = {}
        self.load_g = True
        self.on_g = on_g
        
        # Ensure backward compatibility with older pickled objects
        self.__ensure_required_attributes__()
        
        self.__init_reprdict__()
        self.__init_sgraph__()
        self.__init_paths__(path_data=path_data, path_plot=path_plot)
    
    def __ensure_required_attributes__(self) -> None:
        """
        Ensure all required attributes exist for backward compatibility
        with older pickled objects that may be missing some attributes.
        """
        # Initialize missing core attributes
        required_attrs = {
            'eset': {},
            'fleset': {},
            'lfeset': {},
            'map_node': {},
            'map_edge': {},
            'graph_clustering_utility': None  # Will be initialized later if needed
        }
        
        for attr_name, default_value in required_attrs.items():
            if not hasattr(self, attr_name):
                setattr(self, attr_name, default_value)
        
    #
    def __init_weights__(
        self, 
        values: Union[float, dict] = 1., 
        on_g: str = None
    ) -> None:
        on_g = on_g or self.on_g
        nx.set_edge_attributes(self.gr[on_g], values, 'weight')
        self.upd_GraphRepr_All(on_g)
    #
    def __init_sgraph__(self, init_weights_val: Union[float, dict] = 1.) -> None:
        on_g = self.on_g
        self.eset[on_g] = set(list(self.gr[on_g].edges()))
        if self.load_g:
            edges_data = self.gr[on_g].edges(data=True)
            self.fleset[on_g] = set([
                (u, v) for u, v, _ in edges_data if _.get('weight', 1) < 0
            ])
            self.pflip = self.Ne_n / self.Ne
            self.Ne_flips = int(self.pflip * self.Ne)
            self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
        else:
            self.nflip = int(self.pflip * self.N) # interessante quantita,
            self.Ne_flips = int(self.pflip * self.Ne)
            # i nodi toccati (come variano rispetto alla statistica, cioe 
            # in media sono p * N, ma le fluttuazioni? sono gaussiane...)
            self.fleset[on_g] = set(
                self.get_random_links(self.Ne_flips, on_g=on_g)
            )
            self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
            self.__init_weights__(init_weights_val)
        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices()
    #
    def _verify_pflip(self, pflip: float) -> None:
        """
        Verify if the provided pflip value is within the valid range.

        Parameters
        ----------
        pflip : float
            The probability of flipping edges in the signed graph.

        Raises
        ------
        ValueError
            If pflip is not within the range [LB_PFLIP, UB_PFLIP].
        """
        if not is_in_range(pflip, LB_PFLIP, UB_PFLIP):
            raise ValueError(SG_ERRMSG_PFLIP)
        else:
            self.pflip = pflip
            from ...config.funcs import peq_fstr
            self.peq_str = peq_fstr(pflip)
    #
    # export graph tools
    #
    from ._exports import __export_graph__
    from ._exports import _export_edgel_bin
    from ._exports import _export_eigV 
    from ._exports import export_eigV_all
    #
    # load graph tools
    #
    from ._loaders import __load_graph__
    from ._loaders import _load_eigV
    from ._loaders import load_eigV_all
    #
    def export_ising_clust(self, NoClust: int = 1, exName: str = '') -> None: 
        self.clPname = []
        for i in range(NoClust):
            exname = '_'.join([self.peq_str, exName]) if exName else self.peq_str
            fname = join_non_empty('_', f"cl{i}", exname)+BIN
            self.clPname.append(self.path_ising / fname)
            with open(self.clPname[i], "wb") as f:
                np.array(list(self.biggestClSet[i])).astype(int).tofile(f)
    #
    def remove_ising_clust_files(self):
        for pname in self.clPname:
            os.remove(pname)
    #
    def remove_edgl_file(self):
        os.remove(self.path_exp_edgl)
    #
    def remove_eigV_file(self):
        os.remove(self.eigVPname)
    #
    def remove_exported_files(self):
        if hasattr(self, "path_exp_edgl"):
            self.remove_edgl_file()
        if hasattr(self, "eigVPname"):
            self.remove_eigV_file()
        if hasattr(self, "clPname"):
            self.remove_ising_clust_files()
    #        

    # 
    def set_edgel_from_bin(self, file_path, mode='numpy', on_g: str = SG_GRAPH_REPR):
        match mode:
            case 'numpy':
                dtype = np.dtype([('i', np.uint64), ('j', np.uint64), ('w_ij', np.float64)])
                edge_array = np.fromfile(file_path, dtype=dtype)
                edges = [(int(edge['i']), int(edge['j']), float(edge['w_ij'])) for edge in edge_array]
            case 'struct':
                edges = []
                with open(file_path, "rb") as f:
                    while chunk := f.read(24):  # 2 * uint64 (8 bytes each) + 1 * float64 (8 bytes) = 24 bytes
                        i, j, w_ij = struct.unpack("QQd", chunk)
                        edges.append((i, j, w_ij))
            case _:
                raise ValueError("Invalid mode. Choose 'numpy' or 'struct'.")
        self.gr[on_g].add_weighted_edges_from(edges)
        self.upd_edge_sets(on_g)
        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices(on_g)
    #
    # graph get attributes
    #
    def get_p_fname(
            self, 
            who: str, 
            out_suffix: str = '', 
            ext: str = BIN, 
    ) -> str:
        """Get the file name for exporting or importing graph data."""
        from ...config.funcs import build_p_fname
        return build_p_fname(who, self.pflip, out_suffix=out_suffix, ext=ext)

    def get_node_attributes(
            self, 
            attr: str = 'pos', 
            on_g: str = SG_GRAPH_REPR
    ):
        return nx.get_node_attributes(self.gr[on_g], attr)
    #
    def get_edge_data(self, u: Any, v: Any, thedata: str = 'weight',
                      on_g: str = SG_GRAPH_REPR):
        edge_data = self.gr[on_g].get_edge_data(u, v)
        if edge_data is None:
            raise KeyError(f"Edge ({u}, {v}) not found in graph {on_g}")
        return edge_data[thedata]
    #
    def get_edge_mapping(self, edge, target_g, on_g: str = SG_GRAPH_REPR):
        try:
            return self.map_edge[target_g][on_g][edge]
        except KeyError:
            rev_edge = edge[::-1]
            return self.map_edge[target_g][on_g].get(rev_edge, None)
    #
    def get_edge_color(self, pec: ColorType = "blue", nec: ColorType = "red",
                    thedata: str = 'weight', on_g: str = SG_GRAPH_REPR,
                    continuous: bool = False, cmap: str = 'coolwarm'):
        def map_values(value):
            if continuous:
                norm_value = (value - min_val) / (max_val - min_val)  # Normalize value to [0, 1]
                return get_cmap(cmap)(norm_value)  # Get color from colormap
            return nec if value == -1 else pec if value == 1 else value
        
        arr = nx.get_edge_attributes(self.gr[on_g], thedata)

        if continuous:
            values = np.array(list(arr.values()))
            min_val, max_val = values.min(), values.max()

        return list(map(map_values, arr.values()))
    #  
    def get_graph_neighbors(self, node: Any, on_g: str = SG_GRAPH_REPR):
        return list(self.gr[on_g].neighbors(node))
    #
    def get_adjacency_matrix(self, on_g: str = SG_GRAPH_REPR, 
                             weight: str = 'weight', format: str = 'csr'):
        return nx.to_scipy_sparse_array(self.gr[on_g], 
                                        weight=weight, format=format)
    #
    def get_degree_matrix(self, format: str = 'csr'):
        return spdiags(self.adj.sum(axis=1), 0, *self.adj.shape, format=format)
    #
    def get_abs_degree_matrix(self, format: str = 'csr'):
        return spdiags(abs(self.adj).sum(axis=1), 0, *self.adj.shape, 
                       format=format)
    #
    def get_laplacian(self):
        return self.degm - self.adj
    #
    def get_signed_laplacian(self):
        return self.sdeg - self.adj
    #
    def get_random_links(self, n: int = 1, only_in: str = '',
                         on_g: str = SG_GRAPH_REPR):
        match only_in:
            case ''|'all':
                return random.sample(tuple(self.eset[on_g]), n)
            case '+'|'positive'|'+1'|'plus':
                return random.sample(tuple(self.lfeset[on_g]), n)
            case '-'|'negative'|'-1'|'minus':
                return random.sample(tuple(self.fleset[on_g]), n)
    #
    def get_eigV(self, which: int = 0, binarize: bool = False):
        if binarize:
            return self.get_eigV_binarized(which)
        else:
            eigV = self.eigV[which].squeeze()
            return flip_to_positive_majority_adapted(eigV).squeeze()
    #
    def get_eigV_check(self, which: int = 0, binarize: bool = False, reshaped: bool = False):
        if not hasattr(self, f"eigV") or which >= len(self.eigV):
            self.compute_k_eigvV(k=which+1)
        if reshaped:
            return self.get_eigV(which, binarize).reshape(*self.syshape)
        else:
            return self.get_eigV(which, binarize)
    #
    def get_eigV_binarized(self, which: int = 0):
        eigV = bin_sign(self.eigV[which].squeeze())
        return flip_to_positive_majority(eigV).squeeze()
    #
    def get_eigV_bin_check(self, which: int = 0, reshaped: bool = False):
        if not hasattr(self, f"eigV") or which >= len(self.eigV):
            self.compute_k_eigvV(k=which+1)
        if reshaped:
            return self.get_eigV_binarized(which).reshape(*self.syshape)
        else:
            return self.get_eigV_binarized(which)
    #
    def get_eigV_bin_check_list(self, custom_slice: slice = slice(None), asarray: bool = False) -> Union[List[int], NDArray]:
        maxidx = (custom_slice.stop - 1) if custom_slice.stop is not None else self.N
        assert maxidx <= self.N, SG_ERRMSG_MAXEIGVIDX
        if not hasattr(self, f"eigV"):
            if custom_slice == slice(None):
                self.compute_laplacian_spectrum_weigV()
            else:
                self.compute_k_eigvV(k=maxidx)
        if asarray:
            return np.array([self.get_eigV_bin_check(i) for i in range(maxidx)])
        return [self.get_eigV_bin_check(_) for _ in range(maxidx)]
    #
    # def get_eigV_NDArray(self, which_slice: slice = slice(None), 
    #                      reshaped: bool = False, binarize: bool = False) -> NDArray:
    #     """
    #     Get a NumPy array of eigenvalues for the signed graph.

    #     Parameters
    #     ----------
    #     which_slice : slice, optional
    #         A slice object to specify which eigenvalues to retrieve (default is all).
    #     reshaped : bool, optional
    #         If True, reshapes the output to the shape of the signed graph (default is False).
    #     binarize : bool, optional
    #         If True, returns binarized eigenvalues (default is False).

    #     Returns
    #     -------
    #     NDArray
    #         A NumPy array containing the requested eigenvalues.
    #     """
    #     eigV_arr = self
    
    def get_sgspect_basis(
            self,
            max_factor: int = 2,
            start: int = 0, 
            step: int = 1,
            normalized: bool = False
    ) -> NDArray:
        """
        Generate a basis for the signed graph spectrum.

        Parameters
        ----------
        max_factor : int, optional
            The maximum factor for generating the basis (default is 2).
        start : int, optional
            The starting index for the basis (default is 0).
        step : int, optional
            The step size for generating the basis (default is 1).
        normalized : bool, optional
            If True, normalizes the eigenvalues (default is False).

        Returns
        -------
        NDArray
            A NumPy array containing the generated basis.
        """
        basis_list = np.array([
            self.get_eigV_check(i, reshaped=False) 
            for i in range(start, int(self.N // max_factor), step)
        ])
        if normalized:
            basis_list = normalize_array(basis_list, axis=1)
        return basis_list
    #
    def get_signed_laplacian_embedding(self, k: int = 2):
        return self.eigV[:k]
    #
    def get_subgraph_from_nodes(self, list_of_nodes, on_g: str = SG_GRAPH_REPR):
        return self.gr[on_g].subgraph(list_of_nodes)
    #
    def get_nodes_subgraph_by_kv(self, k, val, on_g: str = SG_GRAPH_REPR):
        G = self.gr[on_g]
        G_yes, G_no = G.copy(), G.copy()
        predicate = val if callable(val) else lambda x: x == val
        for node, v in G.nodes(data=k):
            if predicate(v):
                G_yes.remove_node(node)
            else:
                G_no.remove_node(node)
        return G_yes, G_no
    #
    def get_bineigV_cluster_sizes(self, which: int = 0, 
                      binarize: bool = True, on_g: str = SG_GRAPH_REPR):
        if not all(f"eigV{which}" in self.gr[on_g].nodes[node] 
                   for node in self.gr[on_g].nodes):
            self.load_eigV_on_graph(which, on_g, binarize)
        if not hasattr(self, "clustersY"):
            self.make_clustersYN(f"eigV{which}", +1, on_g)
        cl_len = sorted(map(len, self.clustersY), reverse=True)
        return cl_len
    #
    def get_cluster_distribution(self, which: int = 0, 
                                  on_g: str = SG_GRAPH_REPR,
                                  binarize: bool = True):
        cl_len = self.get_bineigV_cluster_sizes(which, on_g, binarize)
        dictdist_cluster_sizes = {
            size: cl_len.count(size) for size in set(cl_len)
        }
        return dictdist_cluster_sizes
    #
    def get_ferroAntiferro_regions(self, attr_str: str = 's', 
                                        on_g: str = SG_GRAPH_REPR):
        antiGroup = []
        ferroGroup = []
        graph = self.gr[on_g]
        for node, att in graph.nodes(data=attr_str):
            neighbors = graph.neighbors(node)
            if all(graph.nodes[n][attr_str] == -att for n in neighbors):
                antiGroup.append(node)
            elif all(graph.nodes[n][attr_str] == att for n in neighbors):
                ferroGroup.append(node)
        return ferroGroup, antiGroup
    #
    # set graph properties
    #
    def set_node_attributes(self, values: Any, attribute_name: Any, on_g: str = SG_GRAPH_REPR):
        nx.set_node_attributes(self.gr[on_g], values, attribute_name)
    #
    # update graph methods
    #
    def upd_graph_matrices(self, format: str = 'csr', 
                           on_g: str = SG_GRAPH_REPR):
        self.adjacency_matrix = self.get_adjacency_matrix(on_g=on_g)
        self.degree_matrix = self.get_degree_matrix()
        self.signed_degree_matrix = self.get_abs_degree_matrix()
        self.laplacian_matrix = self.get_laplacian()
        self.signed_laplacian_matrix = self.get_signed_laplacian()
        self.upd_Degree(on_g)
    #
    def upd_edge_sets(self, on_g: str = SG_GRAPH_REPR):
        self.eset[on_g] = set(self.gr[on_g].edges())
        self.fleset[on_g] = set([
            (u, v) 
            for u, v, data in self.gr[on_g].edges(data=True)
            if data.get('weight', 1) < 0
        ])
        self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
    #
    def upd_Degree(self, on_g: str = SG_GRAPH_REPR):
        self.degrees = list(dict(self.gr[on_g].degree).values())
    #
    def zip_reprNodes(self, x, on_g: str = SG_GRAPH_REPR):
        return dict(zip(self.gr[on_g], self.gr[x]))
    #
    def zip_reprEdges(self, x, on_g: str = SG_GRAPH_REPR):
        return dict(zip(self.gr[on_g].edges(), self.gr[x].edges()))
    #
    def upd_NodeMap(self, on_g: str = SG_GRAPH_REPR):
        self.map_node[on_g] = {x: 
            {v: k for k, v in self.zip_reprNodes(x, on_g).items()} 
            for x in self.graph_reprs if x != on_g}
    #
    def upd_EdgeMap(self, on_g: str = SG_GRAPH_REPR):
        graph = self.gr[on_g]
        egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
        self.map_edge[on_g] = {x: {v: k for k, v in egraph_to(x).items()} 
                                   for x in self.graph_reprs if x != on_g}
    #
    def upd_GraphRelabel(
            self, 
            to_graph: str = SG_GRAPHINT_REPR, 
            on_g: str = SG_GRAPH_REPR
    ) -> None:
        node_map = self.map_node[to_graph][on_g]
        self.gr[to_graph] = nx.relabel_nodes(self.gr[on_g], node_map)
        self.eset[to_graph] = {x for e in self.eset[on_g] 
                            if (x := self.get_edge_mapping(e, to_graph, on_g))}
        self.fleset[to_graph] = {x for e in self.fleset[on_g]
                            if (x := self.get_edge_mapping(e, to_graph, on_g))}
        self.lfeset[to_graph] = {x for e in self.lfeset[on_g]
                            if (x := self.get_edge_mapping(e, to_graph, on_g))}
        setattr(self, to_graph, self.gr[to_graph])
    #
    def upd_ReprMaps(self, on_g: str = SG_GRAPH_REPR):
        graph = self.gr[on_g]
        ngraph_to = lambda x: dict(zip(graph, self.gr[x]))
        self.map_node[on_g] = {x: {v: k for k, v in ngraph_to(x).items()} 
                                   for x in self.graph_reprs if x != on_g}
        egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
        self.map_edge[on_g] = {x: {v: k for k, v in egraph_to(x).items()} 
                                   for x in self.graph_reprs if x != on_g}
    #
    def upd_GraphRepr_All(self, on_g: str = SG_GRAPH_REPR, 
                          also_itself: bool = True):
        if also_itself:
            self.upd_ReprMaps(on_g=on_g)
        for i in self.graph_reprs:
            if i != on_g:
                self.upd_ReprMaps(on_g=i)
                self.upd_GraphRelabel(to_graph=i, on_g=on_g)
    #
    # graph operations
    #
    def check_Ne_flips(self):
        if self.Ne_flips < 1: raise NflipError(SG_ERRMSG_NFLIP)
    #
    def flip_sel_edges(
            self, 
            links: Any, 
            name: str = 'weight', 
            on_g: str = SG_GRAPH_REPR
    ) -> None:
        neg_weights_dict = {
            (u, v): -1 * self.get_edge_data(u, v, on_g=on_g) 
            for u, v in links
        }
        nx.set_edge_attributes(
            self.gr[on_g],
            values=neg_weights_dict, 
            name=name
        )

        self.fleset[on_g].update(links)
        self.lfeset[on_g].difference_update(links)

        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices(on_g)
    #
    def flip_random_fract_edges(
        self, 
        pflip: float = None, 
        on_g: str = SG_GRAPH_REPR
    ) -> None:
        try:
            if pflip:
                self.pflip = pflip
                self.Ne_flips = int(self.pflip * self.Ne)
                self.check_Ne_flips()
                self.flip_sel_edges(
                    self.get_random_links(self.Ne_flips, on_g), 
                    on_g=on_g
                )
            else:
                self.Ne_flips = int(self.pflip * self.Ne)
                self.check_Ne_flips()
                self.flip_sel_edges(self.fleset[on_g], on_g=on_g)
        except NflipError:
            logger.error(SG_ERRMSG_NFLIP)
            pass
    #
    def unflip_all(self, on_g: str = SG_GRAPH_REPR):
        self.flip_sel_edges(1, on_g=on_g)
    #
    def set_edges_random_normal(self, mu: float = 1.0, sigma: float = 1.0, 
                                 on_g: str = SG_GRAPH_REPR):
        weights = {
            edge: random.normalvariate(mu, sigma) 
            for edge in self.gr[on_g].edges()
        }
        self.set_edge_weights_wij(weights, on_g)
    
    #
    def load_vec_on_nodes(self, vec: NDArray, attr: str,
                          on_g: str = SG_GRAPH_REPR):
        vecNodeAttr = {nd: v for v, nd in zip(vec, self.gr[on_g].nodes)}
        nx.set_node_attributes(self.gr[on_g], values=vecNodeAttr, name=attr)
    #
    def load_eigV_on_graph(self, which: int = 0, on_g: str = SG_GRAPH_REPR, 
                           binarize: bool = False):
        if binarize: 
            eigV = self.get_eigV_bin_check(which=which)
        else:
            try:
                eigV = self.eigV[which]
            except (IndexError, AttributeError):
                self.compute_k_eigvV(k=which+1)
                eigV = self.eigV[which]
        self.load_vec_on_nodes(eigV, f"eigV{which}", on_g)
        # eigV_val_nd = {nd: v for v, nd in zip(eigV, self.gr[on_g].nodes)}
        # nx.set_node_attributes(self.gr[on_g], eigV_val_nd, f"eigV{which}")
    #
    # computations
    #
    def compute_laplacian_spectrum(self, typf: type = np.float64):
        # NOT WORKING, NEEDS FIX
        self.eigv = np.linalg.eigvalsh(self.laplacian_matrix.astype(typf).toarray())
    #
    def compute_laplacian_spectrum_weigV(self, typf: type = np.float64, 
                                         transpose: bool = True, 
                                         with_routine: str = 'numpy'):
        logger.info("Computing eigenvectors for the signed graph Laplacian.")
        match with_routine:                
            case 'cupy':
                slp = cp.asarray(self.slp.astype(typf).toarray())
                self.eigv, self.eigV = cp.linalg.eigh(slp)
                self.eigv, self.eigV = self.eigv.get(), self.eigV.get()
            case 'numpy' | _:
                slp = self.slp.astype(typf).toarray()
                self.eigv, self.eigV = np.linalg.eigh(slp)
        if transpose:
            self.make_eigV_transposed()
    #
    def compute_k_eigvV(self, k: int = 1, with_routine: str = 'scipy',
        which: str = 'SM', typf: type = np.float64, transpose: bool = True):
        if (with_routine in ['numpy', 'cupy']) or k > self.N//2:
            self.compute_laplacian_spectrum_weigV(typf, transpose, with_routine)
        elif with_routine.startswith('scipy'):
            mode = with_routine.split('_')
            mode = mode[-1] if len(mode) > 1 else 'caley'
            self.eigv, self.eigV = scsp_eigsh(
                self.slp.astype(typf), k=k, which=which, mode=mode
            )

            self.make_eigV_transposed()
    #
    def compute_pinf(self, which: int = 0, on_g: str = SG_GRAPH_REPR):
        clustd = np.array(self.get_bineigV_cluster_sizes(which, on_g))
        mclust = clustd[0]
        self.Pinf = mclust / self.N
        self.Pinf_var = np.sum(clustd@clustd-mclust**2)/(np.sum(clustd)-mclust)
        if hasattr(self, "Pinf_dict"):
            self.Pinf_dict[which] = (self.Pinf, self.Pinf_var)
        else:
            self.Pinf_dict = {which: (self.Pinf, self.Pinf_var)}
    #
    def compute_rbim_energy_eigV(self, which: int = 0, use_gpu: bool = False, on_g: str = SG_GRAPH_REPR):
        spins = self.get_eigV_bin_check(which)
        edges = list(self.gr[on_g].edges(data='weight'))
        if not hasattr(self, "energy_eigV_RBIM"):
            self.energy_eigV_RBIM = {}
        self.energy_eigV_RBIM[which] = compute_ising_pairwise_energy(spins, edges, use_gpu=use_gpu)
    #
    def compute_rbim_energy_eigV_all(self, on_g: str = SG_GRAPH_REPR,
                                     **kw_laplspect):
        if not hasattr(self, "eigV"):
            self.compute_laplacian_spectrum_weigV(**kw_laplspect)
        if not hasattr(self, "energy_eigV_RBIM"):
            self.energy_eigV_RBIM = {}
        for which in range(len(self.eigV)):
            if which not in self.energy_eigV_RBIM:
                self.compute_rbim_energy_eigV(which, on_g)
    #
    def get_rbim_energy_eigV(self, which: int = 0):
        if not hasattr(self, "energy_eigV_RBIM"):
            self.compute_rbim_energy_eigV(which)
        elif which not in self.energy_eigV_RBIM.keys():
                self.compute_rbim_energy_eigV(which)
        return self.energy_eigV_RBIM[which]
    #
    def get_all_rbim_energy_eigV(self, as_dict: bool = False, on_g: str = SG_GRAPH_REPR,
                                 **kw_laplspect):
        if not hasattr(self, "energy_eigV_RBIM"):
            self.compute_rbim_energy_eigV_all(on_g, **kw_laplspect)
        if as_dict:
            return self.energy_eigV_RBIM
        else:
            return np.array(list(self.energy_eigV_RBIM.values()))
    #
    # make methods
    #
    def make_eigV_transposed(self):
        self.eigV = self.eigV.T
    #
    def make_rescaled_signed_laplacian(self, MODE: str = 'field'):
        if MODE == 'field':
            self.resLp = self.slp - self.eigv[0] * scsp_identity(self.N)
        elif MODE == 'double':
            from scipy.linalg import eigvalsh

            self.resLp = self.slp - np.array([self.eigv[0]])
            new_eigv0 = eigvalsh(
                self.resLp.astype(np.float64), subset_by_index=[0, 0]
            )
            self.resLp = self.resLp - new_eigv0 * np.identity(self.N)
    #
    def make_graphYN(self, k, val: ConditionalPartitioning, on_g: str = SG_GRAPH_REPR):
        self.graph_clustering_utility[k][val.key][on_g] = self.get_nodes_subgraph_by_kv(k, val.cond_func, on_g)
    #
    def make_clustersYN(self, k, val: ConditionalPartitioning, on_g: str = SG_GRAPH_REPR):
        try:
            graphY, graphN = self.graph_clustering_utility[k][val.key][on_g]
        except:
            self.make_graphYN(k, val, on_g)
            graphY, graphN = self.graph_clustering_utility[k][val.key][on_g]
        #
        self.clustersY = list(nx.connected_components(graphY))
        self.clustersN = list(nx.connected_components(graphN))
        #
        self.numClustersY = len(self.clustersY)
        self.numClustersN = len(self.clustersN)
        #
        if not self.clustersY:  
            warnings.warn(SG_WARNMSG_NOCLUST, SignedGraphWarning)
            self.biggestClSet = []
            self.numClustersBig = 0
            self.gc = None  # Set to None or a default value
            return
        #
        self.biggestClSet = sorted(self.clustersY, key=len, reverse=True)
        self.numClustersBig = len(self.biggestClSet)

        self.gc = max(self.biggestClSet, key=len)
    #
    #
    def make_eigVclustersYN(self, val: ConditionalPartitioning, which: int = 0, 
                            on_g: str = SG_GRAPH_REPR, binarize: bool = True):
        self.load_eigV_on_graph(which, on_g, binarize)
        self.make_clustersYN(f'eigV{which}', val, on_g)
    #
    def make_connected_component_by_edge(self, edge_attr='weight', value=1, 
                                         on_g: str = SG_GRAPH_REPR):
        connected_components = nx.connected_components
        value_edges = [
            (u, v) for u, v, attrs in self.gr[on_g].edges(data=True)
            if attrs.get(edge_attr) == value
        ]
        if not value_edges:
            raise ValueError("No positive edges found in the graph.")
        G_pos = self.gr[on_g].edge_subgraph(value_edges).copy()
        all_connected_components = list(connected_components(G_pos))
        if not all_connected_components:
            raise ValueError("No connected components found with positive edges.")
        self.largest_cc = max(all_connected_components, key=len)
        self.largest_cc_subgraph = G_pos.subgraph(self.largest_cc).copy()
    #
    # cleaners
    #
    def clean_gclutil(self, k, val, on_g: str = SG_GRAPH_REPR):
        self.graph_clustering_utility[k][val][on_g] = None
    #
    # handlers
    #
    def handle_no_clust(self, NoClust: int) -> int | None:
        try:
            if self.numClustersBig == 0:
                logger.error(SG_ERRMSG_ZEROCLUST)
                raise NoClustError(SG_ERRMSG_ZEROCLUST)
            elif NoClust > self.numClustersBig:
                logger.error(SG_ERRMSG_NOCLUST_GTR_NCB)
                raise NoClustError(SG_ERRMSG_NOCLUST_GTR_NCB)
        except NoClustError as exception:
            logger.error(str(exception))
            if str(exception) == SG_ERRMSG_NOCLUST_GTR_NCB:
                logger.error(f"Returning the number of biggest clusters: {len(self.biggestClSet)}")
                return len(self.biggestClSet)  # Correct NoClust
            elif str(exception) == SG_ERRMSG_ZEROCLUST:
                logger.error("No clusters found, returning None.")
                return None  # Indicate the caller should skip
        return NoClust  # If no errors, return the original value