import networkx as nx
import numpy as np
import pickle as pk
#
from functools import partial
from typing import Callable, Any
#
from ...config.const import *
from ...config.funcs import build_p_fname

#
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
        The name of the file to import the eigenvalues from. If not provided,
        it defaults to an empty string.
    ext : str, optional
        The format of the file to import the eigenvalues from. Defaults to '.bin'.
    binarize : bool, optional
        If True, the eigenvalues will be binarized before being stored. Defaults to True.
                
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
