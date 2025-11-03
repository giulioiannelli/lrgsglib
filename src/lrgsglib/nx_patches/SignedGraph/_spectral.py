import logging
import networkx as nx
import numpy as np
from typing import List, Union, Optional, TYPE_CHECKING
from numpy.typing import NDArray
from scipy.sparse import identity as scsp_identity
from scipy.sparse.linalg import eigsh as scsp_eigsh

from ...config.const import *
from ...config.errwar import SignedGraphWarning
from ...utils.basic import (
    flip_to_positive_majority_adapted,
    bin_sign,
    normalize_array,
    dtype_numerical_precision,
)
# info-theory helpers moved to _infotheory

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


#
# Laplacian getters
#
    

#
# Eigenvector helpers
#
def get_eigV(self: "SignedGraph", which: int = 0, binarize: bool = False):
    if binarize:
        return get_eigV_binarized(self, which)
    else:
        eigV = self.eigV[which].squeeze()
        return flip_to_positive_majority_adapted(eigV).squeeze()

def get_eigV_check(
    self,
    which: int = 0,
    binarize: bool = False,
    reshaped: bool = False,
):
    if not hasattr(self, "eigV") or which >= len(self.eigV):
        compute_k_eigvV(self, k=which + 1)
    if reshaped:
        return get_eigV(self, which, binarize).reshape(*self.syshape)
    else:
        return get_eigV(self, which, binarize)

def get_eigV_binarized(self: "SignedGraph", which: int = 0):
    eigV = bin_sign(self.eigV[which].squeeze())
    return flip_to_positive_majority_adapted(eigV).squeeze()

def get_eigV_bin_check(self: "SignedGraph", which: int = 0, reshaped: bool = False):
    if not hasattr(self, "eigV") or which >= len(self.eigV):
        compute_k_eigvV(self, k=which + 1)
    if reshaped:
        return get_eigV_binarized(self, which).reshape(*self.syshape)
    else:
        return get_eigV_binarized(self, which)

def get_eigV_bin_check_list(
    self, custom_slice: slice = slice(None), asarray: bool = False
) -> Union[List[int], NDArray]:
    maxidx_raw = custom_slice.stop if custom_slice.stop else self.N
    maxidx = maxidx_raw - 1 if custom_slice.stop else self.N
    assert maxidx <= self.N, SG_ERRMSG_MAXEIGVIDX
    if not hasattr(self, "eigV"):
        if custom_slice == slice(None):
            compute_laplacian_spectrum_weigV(self)
        else:
            compute_k_eigvV(self, k=maxidx)
    if asarray:
        return np.array([get_eigV_bin_check(self, i) for i in range(maxidx)])
    return [get_eigV_bin_check(self, _) for _ in range(maxidx)]

def get_sgspect_basis(
    self,
    max_factor: int = 2,
    start: int = 0,
    step: int = 1,
    normalized: bool = False,
) -> NDArray:
    basis_list = np.array(
        [get_eigV_check(self, i, reshaped=False) for i in range(start, int(self.N // max_factor), step)]
    )
    if normalized:
        basis_list = normalize_array(basis_list, axis=1)
    return basis_list

def load_vec_on_nodes(self: "SignedGraph", vec: NDArray, attr: str, on_g: str = SG_REPR):
    vecNodeAttr = {nd: v for v, nd in zip(vec, self.gr[on_g].nodes)}
    nx.set_node_attributes(self.gr[on_g], values=vecNodeAttr, name=attr)

def load_eigV_on_graph(self: "SignedGraph", which: int = 0, on_g: str = SG_REPR, binarize: bool = False):
    if binarize:
        eigV = get_eigV_bin_check(self, which=which)
    else:
        try:
            eigV = self.eigV[which]
        except (IndexError, AttributeError):
            compute_k_eigvV(self, k=which + 1)
            eigV = self.eigV[which]
    load_vec_on_nodes(self, eigV, f"eigV{which}", on_g)

#
# Spectral computations
#
def compute_laplacian_spectrum(self: "SignedGraph", typf: type = np.float64):
    # NOT WORKING, NEEDS FIX
    lapl_arr = self.laplacian_matrix.astype(typf).toarray()
    self.eigv = np.linalg.eigvalsh(lapl_arr)

def compute_laplacian_spectrum_weigV(
    self, typf: type = np.float64, transpose: bool = True, backend: str = 'numpy'
):
    logger.info("Computing eigenvectors for the signed graph Laplacian.")
    cached_eigv = getattr(self, "eigv", None)
    cached_eigV = getattr(self, "eigV", None)
    cached_ready = (
        cached_eigv is not None
        and cached_eigV is not None
        and len(cached_eigv) == self.N
        and cached_eigV.shape == (self.N, self.N)
    )
    if cached_ready:
        logger.debug("Using cached signed Laplacian spectrum.")
        eigv_transposed = getattr(self, "_eigV_is_transposed", False)
        if transpose and not eigv_transposed:
            make_eigV_transposed(self)
        elif not transpose and eigv_transposed:
            make_eigV_column_major(self)
        return

    match backend:
        case 'cupy':
            import cupy as cp
            slp = cp.asarray(self.slp.astype(typf).toarray())
            self.eigv, self.eigV = cp.linalg.eigh(slp)
            self.eigv, self.eigV = self.eigv.get(), self.eigV.get()
        case 'numpy' | _:
            slp = self.slp.astype(typf).toarray()
            self.eigv, self.eigV = np.linalg.eigh(slp)

    self._eigV_is_transposed = False
    if transpose:
        make_eigV_transposed(self)

def compute_k_eigvV(
    self, k: int = 1, backend: str = 'scipy', which: str = 'SM', typf: type = np.float64, transpose: bool = True
):
    if (backend in ['numpy', 'cupy']) or k > self.N // 2:
        compute_laplacian_spectrum_weigV(self, typf, transpose, backend)
    elif backend.startswith('scipy'):
        mode = backend.split('_')
        mode = mode[-1] if len(mode) > 1 else 'caley'
        self.eigv, self.eigV = scsp_eigsh(self.slp.astype(typf), k=k, which=which, mode=mode)

        self._eigV_is_transposed = False
        make_eigV_transposed(self)

    

    

#
# Layout flip helpers
#
def make_eigV_transposed(self: "SignedGraph"):
    if getattr(self, "_eigV_is_transposed", False):
        return
    self.eigV = self.eigV.T
    self._eigV_is_transposed = True

def make_eigV_column_major(self: "SignedGraph"):
    if not getattr(self, "_eigV_is_transposed", False):
        return
    self.eigV = self.eigV.T
    self._eigV_is_transposed = False

    
