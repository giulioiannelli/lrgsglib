import struct
#
import numpy as np
import pickle as pk
import networkx as nx
#
from typing import Union, Callable, Any
from pathlib import Path
from functools import partial
#
from ...config.const import *
from ...config.funcs import build_p_fname
from ...utils.basic import join_non_empty
#
def __export_graph__(self, file_name: str = '', export_mode: str = SG_EXPORT_M):
    match export_mode:
        case 'pkl'|'pk'|'pickle':
            fname = file_name or self.std_fname + PKL 
            self.path_graph_graph = self.path_graph / fname
            pk.dump(self.G, open(self.path_graph_graph, "wb"), 
                    pk.HIGHEST_PROTOCOL)
        case 'gml':
            fname = self.std_fname + GML
            self.path_graph_graph = self.path_graph / fname
            nx.write_gml(self.G, self.path_graph_graph)
        case 'graphml':
            fname = self.std_fname + XML
            self.path_graph_graph = self.path_graph / fname
            nx.write_graphml(self.G, self.path_graph_graph)
#
def __export_data_tofile__(
        self, 
        who: str,
        export_func: Callable,
        path_sgdata: Path,
        file_name: str = '',
        exName: str = '', 
        ext: str = '', 
        **kwargs: Any
) -> None:
    """Load data from a file into the SignedGraph instance."""
    if not hasattr(self, who):
        raise AttributeError(f"SignedGraph instance has no attribute '{who}'")
    fname = file_name or self.get_p_fname(who=who, out_suffix=exName, ext=ext)
    export_func(self, path_sgdata / fname, **kwargs)
#
def _export_eigV(
        self,
        path: Path,
        binarize: bool = True,
        ext: str = '.bin'
) -> None:
    """Export eigenvalues to a file."""
    dtype = np.int8 if binarize else np.float64
    if not hasattr(self, "eigV"):
        self.compute_laplacian_spectrum_weigV()
    if binarize:
        outarr = self.get_eigV_bin_check_list(asarray=True).astype(dtype)
    else:
        outarr = self.eigV.astype(dtype)
    match ext:
        case '.bin':
            with open(path, "wb") as f:
                outarr.tofile(f)
        case '.txt':
            with open(path, "w") as f:
                for i in range(outarr.shape[0]):
                    np.savetxt(f, outarr[i], fmt='%.3g')
        case '.npz':
            np.savez(path, eigV=outarr)
        case '.pkl'|'pickle':
            with open(path, "wb") as f:
                pk.dump(outarr, f, pk.HIGHEST_PROTOCOL)
        case _:
            raise ValueError(f"Unsupported format: {ext}. Supported formats\
                              are: .bin, .txt, .npz, .pkl")
#
def export_eigV_all(
    self, 
    exName: str = '', 
    ext: str = '.bin', 
    binarize: bool = True, 
    path_sgdata: Path = None
) -> None:
    __export_data_tofile__(
        self,
        who='eigV',
        export_func=partial(_export_eigV, binarize=binarize),
        path_sgdata=path_sgdata or self.path_graph,
        exName=exName,
        ext=ext,
    )

def _export_edgel_bin(
        self,
        exName: str = '',
        mode: str = 'numpy',
        on_g: str = SG_GRAPH_REPR
) -> None:
    fname = build_p_fname('edgelist', self.pflip, out_suffix=exName, ext=BIN)
    self.path_exp_edgl = self.path_graph / fname
    #
    edges = self.gr[on_g].edges(data='weight')
    match mode:
        case 'numpy':
            dtype = [('i', np.uint64), ('j', np.uint64), ('w_ij', np.float64)]
            edge_array = np.array(list(edges), dtype=dtype)
            with open(self.path_exp_edgl, "wb") as f:
                edge_array.tofile(f)
        case 'struct':
            with open(self.path_exp_edgl, "wb") as f:
                for edge in edges:
                    assert len(edge) == 3, "Edge must be: (i, j, w_ij)"
                    f.write(struct.pack("QQd", *edge))
#
# def export_adj_bin(self, verbose: bool = False) -> None:
#     rowarr = [row[i:] for i, row in enumerate(self.adjacency_matrix.toarray())]
#     exname = self.path_graph / f"adj_{self.std_fname}.bin"
#     self.adjFile = open(exname, "wb")
#     with self.adjFile as f:
#         for i in range(len(rowarr)):
#             rowarr[i].astype("float64").tofile(f)
#
# def _export_eigV_core(self, f: Union[str, Path], which: int = 0, binarize: bool = True) -> None:
#     if binarize:
#         self.get_eigV_bin_check(which).astype("int8").tofile(f)
#     else:
#         self.eigV[which].astype("float64").tofile(f)
#
# def _export_eigV(self, which: int = 0,  exName: str = '', verbose: bool = False, binarize: bool = True) -> None:
#     # add support for different formats
#     fname =  build_p_fname(f'eigV{which}', self.pflip, out_suffix=exName, ext=BIN)
#     self.eigVPname = self.path_graph / fname
#     #
#     with open(self.eigVPname, "wb") as f:
#         _export_eigV_core(self, f, which, binarize)
#
# def _export_eigV_all(self, exName: str = '', ext: str = '.bin', binarize: bool = True, custom_path: Path = '') -> None:
#     if not hasattr(self, "eigV"):
#         self.compute_laplacian_spectrum_weigV()
#     #
#     if custom_path:
#         self.path_exp_eigV = custom_path
#     else:
#         fname = build_p_fname('eigV', self.pflip, out_suffix=exName, ext=ext)
#         self.path_exp_eigV = self.path_graph / fname
#     #
#     if binarize:
#         outarr = self.get_eigV_bin_check_list(asarray=True).astype("int8")
#     else:
#         outarr = self.eigV.astype("float64")


