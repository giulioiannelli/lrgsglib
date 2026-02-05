import networkx as nx
import numpy as np
import pickle as pk
import struct
#
from functools import partial
from typing import Callable, Any, TYPE_CHECKING
#
from ....config.const import *
from ....config.funcs import build_p_fname

#
if TYPE_CHECKING:
    from .SignedGraphNX import SignedGraphNX


def __load_graph__(
        self,
        file_name: str = '',
        import_mode: str = SG_LOAD_M
) -> None:
    match import_mode:
        case 'pkl'|'pk'|'pickle':
            fname = file_name or self.std_fname + PKL 
            self.path_load = self.path_graph / fname
            load = pk.load(open(self.path_load, "rb"))
        case 'gml':
            fname = self.std_fname + GML
            self.path_load = self.path_graph / fname
            load = nx.read_gml(self.path_load)
        case 'graphml':
            fname = self.std_fname + XML
            self.path_load = self.path_graph / fname
            nx.read_graphml(self.path_load)
    return load
#
def __load_data_fromfile__(
        self, 
        who: str,
        loader_func: Callable,
        path_sgdata: Path,
        file_name: str = '',
        exName: str = '', 
        ext: str = '', 
        **kwargs: Any
) -> None:
    """Load data from a file into the SignedGraph instance."""
    fname = file_name or self.get_p_fname(who=who, out_suffix=exName, ext=ext)
    load = loader_func(self, path_sgdata / fname, **kwargs)
    setattr(self, f"{who}", load)
#
def _load_eigV(
        self,
        path: Path,
        binarize: bool = True,
        ext: str = '.bin'
) -> None:
    """Load eigenvalues from a file and store them in the eigV attribute."""
    dtype = np.int8 if binarize else np.float64
    if path.exists():
        match ext:
            case '.bin':
                load = np.fromfile(path, dtype=dtype)
            case '.txt':
                load = np.loadtxt(path, dtype=dtype)
            case '.npz':
                load = np.load(path).astype(dtype)
            case '.pkl'|'pickle':
                load = pk.load(open(path, "rb"))
            case _: 
                raise ValueError(f"Unsupported file extension: {ext}")
    else:
        raise FileNotFoundError(f"Eigenvalue file {path} not found.")
    return load
#
def load_eigV_all(
    self, 
    exName: str = '', 
    ext: str = '.bin',
    path_sgdata: Path = None,
    binarize: bool = True,
    reshape: bool = True
) -> None:
    """
    Import eigenvalues from a file and store them in the eigV attribute.
    If the file does not exist, raises a FileNotFoundError.
    
    Parameters
    ----------
    exName : str, optional
        The name of the file to import the eigenvalues from. 
        If not provided, it defaults to an empty string.
    ext : str, optional
        The format of the file to import the eigenvalues from. 
        Defaults to '.bin'.
    binarize : bool, optional
        If True, the eigenvalues will be binarized before being stored. 
        Defaults to True.
                
    Raises
    ------
    FileNotFoundError
        If the eigenvalue file does not exist.
    """
    __load_data_fromfile__(
        self,
        who='eigV', 
        loader_func=partial(_load_eigV, binarize=binarize),
        path_sgdata=path_sgdata or self.path_graph, 
        exName=exName, 
        ext=ext,
    )
    if reshape:
        self.eigV = self.eigV.reshape(-1, self.N)
#
def set_edgel_from_bin(
        self, 
        file_path, 
        mode='numpy', 
        on_g: str = SG_REPR
):
    """
    Load edge list from binary file and add to graph.
    
    Parameters
    ----------
    file_path : str or Path
        Path to the binary file containing edge data.
    mode : str, default 'numpy'
        Loading mode: 'numpy' or 'struct'.
    on_g : str, default SG_REPR
        Graph representation to update.
        
    Notes
    -----
    The 'numpy' mode expects a structured array with fields 'i', 'j', 
    and 'w_ij'.
    The 'struct' mode reads raw binary data (uint64, uint64, float64).
    """
    match mode:
        case 'numpy':
            dtype = np.dtype([
                ('i', np.uint64), 
                ('j', np.uint64), 
                ('w_ij', np.float64)
            ])
            edge_array = np.fromfile(file_path, dtype=dtype)
            edges = [
                (int(edge['i']), int(edge['j']), float(edge['w_ij'])) 
                for edge in edge_array
            ]
        case 'struct':
            edges = []
            with open(file_path, "rb") as f:
                # 2 * uint64 (8 bytes each) + 1 * float64 (8 bytes)
                # = 24 bytes total
                while chunk := f.read(24):
                    i, j, w_ij = struct.unpack("QQd", chunk)
                    edges.append((i, j, w_ij))
        case _:
            raise ValueError("Invalid mode. Choose 'numpy' or 'struct'.")
    self.gr[on_g].add_weighted_edges_from(edges)
    self.upd_edge_sets(on_g)
    self.upd_GraphRepr_All(on_g)
    self.upd_graph_matrices(on_g)
#
