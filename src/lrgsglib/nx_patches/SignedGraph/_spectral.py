import logging
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
# computation methods
#
def compute_laplacian_spectrum(self: "SignedGraph", typf: type = np.float64):
    """
    Compute eigenvalues only of the signed Laplacian (no eigenvectors).
    
    This function computes only the eigenvalues of the signed Laplacian matrix
    without computing eigenvectors, which is faster and uses less memory when
    eigenvectors are not needed. Uses dense eigenvalue decomposition.
    
    Parameters
    ----------
    typf : type, default np.float64
        Data type for computation (e.g., np.float32 for memory savings,
        np.float64 for higher precision).
    
    Returns
    -------
    None
        Results are stored in instance attribute:
        - self.eigv : ndarray of shape (N,)
            Eigenvalues in ascending order.
    
    Notes
    -----
    This function is useful when you only need spectral properties (like
    spectral gap, trace, determinant) without needing the eigenvectors
    themselves. It's significantly faster and more memory-efficient than
    `compute_laplacian_spectrum_weigV()`.
    
    **Memory Usage:**
    This function only stores N eigenvalues (8N bytes for float64), while
    the full eigendecomposition stores N² eigenvector components plus N
    eigenvalues (8N² + 8N bytes for float64).
    
    **Use Cases:**
    - Computing spectral invariants (trace, determinant, spectral radius)
    - Analyzing spectral gap
    - Studying energy landscapes
    - Memory-constrained environments
    
    See Also
    --------
    compute_laplacian_spectrum_weigV : Compute both eigenvalues and eigenvectors
    compute_k_eigvV : Compute k smallest eigenvalues and eigenvectors
    
    Examples
    --------
    >>> # Compute eigenvalues only
    >>> compute_laplacian_spectrum(graph)
    >>> graph.eigv.shape
    (100,)
    
    >>> # Access spectral properties
    >>> spectral_gap = graph.eigv[1] - graph.eigv[0]
    >>> trace = graph.eigv.sum()
    
    >>> # Use lower precision
    >>> compute_laplacian_spectrum(graph, typf=np.float32)
    """
    lapl_arr = self.slp.astype(typf).todense()
    self.eigv = np.linalg.eigvalsh(lapl_arr)
#
def compute_k_eigvV(
    self: "SignedGraph",
    k: int = 1,
    backend: str = 'scipy',
    transpose: bool = True,
    flip_to_pos: bool = True,
    typf: type = np.float64,
    which: str = 'SM'
):
    """
    Compute k smallest eigenvectors of the signed Laplacian.
    
    Parameters
    ----------
    k : int, default 1
        Number of eigenvectors to compute.
    backend : str, default 'scipy'
        Backend to use: 'scipy', 'numpy', or 'cupy'.
    transpose : bool, default True
        If True, eigenvectors stored as row-major (M×N) where M is number
        of eigenvectors and N is vector length. Access as self.eigV[i].
        If False, stored as column-major (N×M).
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    typf : type, default np.float64
        Data type for computation.
    which : str, default 'SM'
        Which eigenvalues to find ('SM' = smallest magnitude).
    """
    if (backend in ['numpy', 'cupy']) or k > self.N // 2:
        compute_laplacian_spectrum_weigV(
            self, backend, transpose, flip_to_pos, typf
        )
    elif backend.startswith('scipy'):
        mode = backend.split('_')
        mode = mode[-1] if len(mode) > 1 else 'caley'
        self.eigv, self.eigV = scsp_eigsh(
            self.slp.astype(typf), k=k, which=which, mode=mode
        )

        # At this point eigV is column-major (N, k) from scipy
        # Apply flip to positive majority if requested (on columns)
        if flip_to_pos:
            # Use axis=0 to flip each column independently
            self.eigV = flip_to_positive_majority_adapted(self.eigV, axis=0)

        # Set transposed state and apply transpose if requested
        self._eigV_is_transposed = False
        if transpose:
            make_eigV_transposed(self)


def compute_k_adj_eigvV(
    self: "SignedGraph",
    k: int = 1,
    backend: str = 'scipy',
    transpose: bool = True,
    flip_to_pos: bool = True,
    typf: type = np.float64,
    which: str = 'LM'
):
    """
    Compute k eigenvectors of the adjacency matrix with sparse/dense backends.

    Parameters
    ----------
    k : int, default 1
        Number of eigenvectors to compute.
    backend : str, default 'scipy'
        Backend to use: 'scipy', 'numpy', or 'cupy'.
    transpose : bool, default True
        If True, eigenvectors stored as row-major (M×N) where M is number of
        eigenvectors and N is vector length. If False, stored as column-major.
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    typf : type, default np.float64
        Data type for computation.
    which : str, default 'LM'
        Which eigenvalues to find ('LM' = largest magnitude by default).
    """
    if getattr(self, "adj", None) is None:
        raise ValueError(
            "Adjacency matrix is not initialized; build it before computing "
            "adjacency eigenpairs."
        )

    if (backend in ['numpy', 'cupy']) or k > self.N // 2:
        compute_adjacency_spectrum_weigV(
            self, backend=backend, transpose=transpose,
            flip_to_pos=flip_to_pos, typf=typf
        )
        return

    if not backend.startswith('scipy'):
        raise ValueError(
            f"Unsupported backend '{backend}' for sparse adjacency eigensolver. "
            "Use a scipy backend or fall back to dense computation."
        )

    mode = backend.split('_')
    mode = mode[-1] if len(mode) > 1 else 'caley'
    self.adj_eigv, self.adj_eigV = scsp_eigsh(
        self.adj.astype(typf), k=k, which=which, mode=mode
    )

    if flip_to_pos:
        self.adj_eigV = flip_to_positive_majority_adapted(self.adj_eigV, axis=0)

    self._adj_eigV_is_transposed = False
    if transpose:
        make_adj_eigV_transposed(self)
#
def compute_laplacian_spectrum_weigV(
    self: "SignedGraph",
    backend: str = 'numpy',
    transpose: bool = True,
    flip_to_pos: bool = True,
    typf: type = np.float64
):
    """
    Compute full eigendecomposition of the signed Laplacian.
    
    This function computes all eigenvalues and eigenvectors of the signed
    Laplacian matrix using dense eigendecomposition. Results are cached and
    reused if already computed with matching dimensions.
    
    Parameters
    ----------
    backend : str, default 'numpy'
        Backend to use for computation:
        - 'numpy': Use numpy.linalg.eigh (CPU, dense)
        - 'cupy': Use cupy.linalg.eigh (GPU, requires cupy installed)
    transpose : bool, default True
        If True, eigenvectors stored as row-major (M×N) where M is number
        of eigenvectors and N is vector length. Access as self.eigV[i].
        If False, stored as column-major (N×M).
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
        This ensures consistent eigenvector orientation across runs.
    typf : type, default np.float64
        Data type for computation (e.g., np.float32 for memory savings,
        np.float64 for higher precision).
    
    Returns
    -------
    None
        Results are stored in instance attributes:
        - self.eigv : ndarray of shape (N,)
            Eigenvalues in ascending order.
        - self.eigV : ndarray of shape (M, N) or (N, M)
            Eigenvectors corresponding to eigenvalues. Shape depends on
            transpose parameter.
        - self._eigV_is_transposed : bool
            Flag indicating current storage layout.
    
    Notes
    -----
    **Caching Behavior:**
    The function checks if a complete spectrum has already been computed
    (self.eigv and self.eigV exist with correct dimensions). If so, it only
    adjusts the transpose state if needed, avoiding redundant computation.
    
    **Storage Layout:**
    - Column-major (transpose=False): eigV has shape (N, N), access i-th
      eigenvector as eigV[:, i]
    - Row-major (transpose=True): eigV has shape (N, N), access i-th
      eigenvector as eigV[i, :]
    
    Row-major is the default as it provides more intuitive indexing and
    better cache locality for subsequent operations on individual eigenvectors.
    
    **Eigenvector Orientation:**
    Both numpy and cupy return eigenvectors in column-major format. The
    flip_to_pos operation is applied along axis=0 (columns) before any
    transposition. This ensures that each eigenvector has more positive
    than negative components, providing consistent orientation.
    
    **Backend Selection:**
    - Use 'numpy' for general purposes (CPU-based, no dependencies)
    - Use 'cupy' for large graphs where GPU acceleration provides significant
      speedup (requires cupy installation and compatible GPU)
    
    See Also
    --------
    compute_k_eigvV : Compute k smallest eigenvectors (efficient for small k)
    make_eigV_transposed : Convert to row-major storage
    make_eigV_column_major : Convert to column-major storage
    
    Examples
    --------
    >>> # Compute full spectrum with default settings
    >>> compute_laplacian_spectrum_weigV(graph)
    >>> graph.eigv.shape
    (100,)
    >>> graph.eigV.shape
    (100, 100)
    
    >>> # Use GPU acceleration
    >>> compute_laplacian_spectrum_weigV(graph, backend='cupy')
    
    >>> # Use lower precision for memory savings
    >>> compute_laplacian_spectrum_weigV(graph, typf=np.float32)
    
    >>> # Store in column-major format
    >>> compute_laplacian_spectrum_weigV(graph, transpose=False)
    >>> v0 = graph.eigV[:, 0]  # First eigenvector
    """
    logger.info("Computing eigenvectors for the signed graph Laplacian.")
    
    cached_eigv = getattr(self, "eigv", None)
    cached_eigV = getattr(self, "eigV", None)
    cached_ready = (
        cached_eigv is not None
        and cached_eigV is not None
        and len(cached_eigv) == self.N
        and cached_eigV.shape[0] == self.N
        and cached_eigV.shape[1] == self.N
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
            try:
                import cupy as cp
                slp = cp.asarray(self.slp.astype(typf).toarray())
                self.eigv, self.eigV = cp.linalg.eigh(slp)
                self.eigv, self.eigV = self.eigv.get(), self.eigV.get()
            except ImportError:
                logger.warning("cupy not available, falling back to numpy backend.")
                slp = self.slp.astype(typf).toarray()
                self.eigv, self.eigV = np.linalg.eigh(slp)
        case 'numpy' | _:
            slp = self.slp.astype(typf).toarray()
            self.eigv, self.eigV = np.linalg.eigh(slp)

    if flip_to_pos:
        self.eigV = flip_to_positive_majority_adapted(self.eigV, axis=0)

    self._eigV_is_transposed = False
    if transpose:
        make_eigV_transposed(self)


def compute_adjacency_spectrum_weigV(
    self: "SignedGraph",
    backend: str = 'numpy',
    transpose: bool = True,
    flip_to_pos: bool = True,
    typf: type = np.float64
):
    """
    Compute the full eigendecomposition of the adjacency matrix.

    Parameters
    ----------
    backend : str, default 'numpy'
        Backend to use for computation: 'numpy' or 'cupy'.
    transpose : bool, default True
        If True, eigenvectors stored as row-major (M×N). If False, column-major.
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    typf : type, default np.float64
        Data type for computation.
    """
    if getattr(self, "adj", None) is None:
        raise ValueError(
            "Adjacency matrix is not initialized; build it before computing "
            "adjacency eigenpairs."
        )

    logger.info("Computing eigenvectors for the adjacency matrix.")

    cached_eigv = getattr(self, "adj_eigv", None)
    cached_eigV = getattr(self, "adj_eigV", None)
    cached_ready = (
        cached_eigv is not None
        and cached_eigV is not None
        and cached_eigV.shape[0] == self.N
        and cached_eigV.shape[1] == self.N
    )

    if cached_ready:
        logger.debug("Using cached adjacency spectrum.")
        eigv_transposed = getattr(self, "_adj_eigV_is_transposed", False)
        if transpose and not eigv_transposed:
            make_adj_eigV_transposed(self)
        elif not transpose and eigv_transposed:
            make_adj_eigV_column_major(self)
        return

    def _dense_adj():
        adj_matrix = self.adj.astype(typf)
        return adj_matrix.toarray() if hasattr(adj_matrix, "toarray") else np.asarray(adj_matrix, dtype=typf)

    match backend:
        case 'cupy':
            try:
                import cupy as cp
                dense_adj = _dense_adj()
                adj_gpu = cp.asarray(dense_adj)
                eigv, eigV = cp.linalg.eigh(adj_gpu)
                self.adj_eigv, self.adj_eigV = eigv.get(), eigV.get()
            except ImportError:
                logger.warning(
                    "cupy not available, falling back to numpy backend for adjacency spectrum."
                )
                dense_adj = _dense_adj()
                self.adj_eigv, self.adj_eigV = np.linalg.eigh(dense_adj)
        case 'numpy' | _:
            dense_adj = _dense_adj()
            self.adj_eigv, self.adj_eigV = np.linalg.eigh(dense_adj)

    if flip_to_pos:
        self.adj_eigV = flip_to_positive_majority_adapted(self.adj_eigV, axis=0)

    self._adj_eigV_is_transposed = False
    if transpose:
        make_adj_eigV_transposed(self)
#
# get methods
#
def get_eigV_squeezed(self: "SignedGraph", which: int = 0):
    """
    Get a specific eigenvector from the cached spectrum, squeezed to 1D.
    
    This function retrieves the eigenvector at the specified index from the
    precomputed eigenvector matrix. The eigenvector is squeezed to remove
    any singleton dimensions, ensuring it's returned as a 1D array.
    
    Parameters
    ----------
    which : int, default 0
        Index of the eigenvector to retrieve. Must be less than the number
        of computed eigenvectors.
    
    Returns
    -------
    ndarray
        The requested eigenvector as a 1D array.
    
    Raises
    ------
    AttributeError
        If eigenvectors have not been computed yet (eigV attribute missing).
    IndexError
        If the requested index is out of range for the computed eigenvectors.
    
    Notes
    -----
    This function does not trigger computation. Use `get_eigV_check()` if you
    want to automatically compute eigenvectors when they're not available.
    
    Examples
    --------
    >>> eigvec = get_eigV_squeezed(graph, which=2)
    >>> eigvec.shape
    (100,)  # 1D array for a graph with 100 nodes
    """
    if hasattr(self, "eigV"):
        try: 
            return self.eigV[which].squeeze()
        except IndexError:
            raise IndexError("Eigenvector index out of range.")
    else:
        raise AttributeError("eigV not found in SignedGraph instance.")
#
def get_eigV_binarized(self: "SignedGraph", which: int = 0):
    """
    Get a binarized eigenvector with values in {-1, +1}.
    
    This function retrieves an eigenvector and applies binary sign conversion,
    mapping all values to either -1 or +1. Zero values are converted to +1.
    
    Parameters
    ----------
    which : int, default 0
        Index of the eigenvector to retrieve and binarize.
    
    Returns
    -------
    ndarray
        The binarized eigenvector as a 1D array with values in {-1, +1}.
    
    Raises
    ------
    AttributeError
        If eigenvectors have not been computed yet (eigV attribute missing).
    IndexError
        If the requested index is out of range for the computed eigenvectors.
    
    Notes
    -----
    This function applies the `bin_sign` transformation which:
    - Converts negative values to -1
    - Converts positive values to +1
    - Converts zero values to +1 (treating zero as positive)
    
    For automatic computation when eigenvectors are missing, use 
    `get_eigV_bin_check()` instead.
    
    See Also
    --------
    get_eigV_squeezed : Get raw eigenvector without binarization
    get_eigV_bin_check : Binarized eigenvector with automatic computation
    bin_sign : The underlying binarization function
    
    Examples
    --------
    >>> eigvec = get_eigV_binarized(graph, which=1)
    >>> set(eigvec)
    {-1, 1}  # Only -1 and +1 values
    """
    return bin_sign(get_eigV_squeezed(self, which))
#
def get_eigV(
    self: "SignedGraph", 
    which: int = 0, 
    binarize: bool = False,
    reshape: bool = False
):
    """
    Get an eigenvector with optional binarization and reshaping.
    
    This is a flexible accessor function that retrieves an eigenvector with
    optional post-processing. It combines the functionality of getting raw or
    binarized eigenvectors and optionally reshaping them to the graph's
    system shape.
    
    Parameters
    ----------
    which : int, default 0
        Index of the eigenvector to retrieve.
    binarize : bool, default False
        If True, return binarized eigenvector with values in {-1, +1}.
        If False, return raw continuous eigenvector values.
    reshape : bool, default False
        If True, reshape the eigenvector to self.syshape (useful for
        lattice graphs). If False, return as 1D array.
    
    Returns
    -------
    ndarray
        The requested eigenvector, either 1D or reshaped to system dimensions.
        Values are either continuous (binarize=False) or in {-1, +1}
        (binarize=True).
    
    Raises
    ------
    AttributeError
        If eigenvectors have not been computed yet (eigV attribute missing).
    IndexError
        If the requested index is out of range for the computed eigenvectors.
    
    Notes
    -----
    This function does not trigger computation. For automatic computation when
    eigenvectors are missing, use `get_eigV_check()` instead.
    
    The reshape option is particularly useful for lattice graphs where nodes
    are arranged in a regular grid and you want to visualize the eigenvector
    as a 2D/3D array matching the lattice structure.
    
    See Also
    --------
    get_eigV_check : Same functionality with automatic computation
    get_eigV_squeezed : Get raw eigenvector without options
    get_eigV_binarized : Get binarized eigenvector
    
    Examples
    --------
    >>> # Get raw eigenvector
    >>> v0 = get_eigV(graph, which=0)
    >>> v0.shape
    (100,)
    
    >>> # Get binarized eigenvector
    >>> v0_bin = get_eigV(graph, which=0, binarize=True)
    >>> set(v0_bin)
    {-1, 1}
    
    >>> # Get eigenvector reshaped to lattice dimensions
    >>> v0_grid = get_eigV(graph, which=0, reshape=True)
    >>> v0_grid.shape
    (10, 10)  # For a 10x10 lattice
    """
    if binarize:
        eigV = get_eigV_binarized(self, which)
    else:
        eigV = get_eigV_squeezed(self, which)
    if reshape:
        return eigV.reshape(*self.syshape)
    else:
        return eigV
#
def get_eigV_check(
    self: "SignedGraph", 
    which: int = 0, 
    binarize: bool = False, 
    reshaped: bool = False,
    backend: str = 'scipy',
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
):
    """
    Get an eigenvector with automatic computation if not available.
    
    This function extends `get_eigV()` by automatically computing eigenvectors
    when they haven't been computed yet or when the requested index is out of
    range. It intelligently chooses between computing k eigenvectors (for small
    k) or the full spectrum (for large k) based on the requested index.
    
    Parameters
    ----------
    which : int, default 0
        Index of the eigenvector to retrieve. Must be non-negative and less
        than the number of nodes N in the graph.
    binarize : bool, default False
        If True, return binarized eigenvector with values in {-1, +1}.
        If False, return raw continuous eigenvector values.
    reshaped : bool, default False
        If True, reshape the eigenvector to self.syshape (useful for
        lattice graphs). If False, return as 1D array.
    backend : str, default 'scipy'
        Backend to use for computation: 'scipy', 'numpy', or 'cupy'.
        For which < N/2, scipy backends use sparse methods; numpy/cupy use dense.
        For which >= N/2, always uses dense computation with numpy or cupy.
    typf : type, default np.float64
        Data type for computation (e.g., np.float32, np.float64).
    transpose : bool, default True
        If True, eigenvectors stored as row-major (M×N) where M is number
        of eigenvectors and N is vector length. Access as self.eigV[i].
        If False, stored as column-major (N×M).
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    
    Returns
    -------
    ndarray
        The requested eigenvector, either 1D or reshaped to system dimensions.
        Values are either continuous (binarize=False) or in {-1, +1}
        (binarize=True).
    
    Notes
    -----
    This function optimizes computation by:
    - Computing only k eigenvectors when which < N/2 using sparse methods
    - Computing the full spectrum when which >= N/2 using dense methods
    
    The automatic computation ensures you can request any eigenvector without
    manually checking if it has been computed, making the API more convenient
    for interactive use.
    
    The backend parameter allows control over the computation method:
    - 'scipy' or 'scipy_caley': Use scipy.sparse.linalg.eigsh (sparse efficient)
    - 'numpy': Use numpy.linalg.eigh (dense computation)
    - 'cupy': Use cupy.linalg.eigh (GPU-accelerated, requires cupy)
    
    See Also
    --------
    get_eigV : Get eigenvector without automatic computation
    get_eigV_bin_check : Convenience wrapper for binarized eigenvectors
    compute_k_eigvV : Compute k smallest eigenvectors
    compute_laplacian_spectrum_weigV : Compute full spectrum
    
    Examples
    --------
    >>> # Request eigenvector 5 - will compute if needed
    >>> v5 = get_eigV_check(graph, which=5)
    >>> v5.shape
    (100,)
    
    >>> # Request binarized and reshaped eigenvector
    >>> v2_bin = get_eigV_check(graph, which=2, binarize=True, reshaped=True)
    >>> set(v2_bin.ravel())
    {-1, 1}
    >>> v2_bin.shape
    (10, 10)  # For a 10x10 lattice
    
    >>> # Request with GPU acceleration
    >>> v0 = get_eigV_check(graph, which=0, backend='cupy')
    
    >>> # Request with specific precision
    >>> v1 = get_eigV_check(graph, which=1, typf=np.float32)
    """
    if not hasattr(self, "eigV") or which >= len(self.eigV):
        if which < self.N // 2:
            compute_k_eigvV(
                self, k=which + 1, backend=backend, typf=typf,
                transpose=transpose, flip_to_pos=flip_to_pos
            )
        else:
            compute_laplacian_spectrum_weigV(
                self, typf=typf, transpose=transpose,
                backend=backend, flip_to_pos=flip_to_pos
            )
    #
    return get_eigV(self, which=which, binarize=binarize, reshape=reshaped)


def get_eigV_bin_check(
    self: "SignedGraph",
    which: int = 0,
    reshaped: bool = False,
    backend: str = 'scipy',
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
):
    """
    Get a binarized eigenvector with automatic computation if not available.
    
    This is a convenience wrapper around `get_eigV_check()` that always returns
    binarized eigenvectors with values in {-1, +1}. It combines automatic
    computation with binarization in a single function call.
    
    Parameters
    ----------
    which : int, default 0
        Index of the eigenvector to retrieve. Must be non-negative and less
        than the number of nodes N in the graph.
    reshaped : bool, default False
        If True, reshape the eigenvector to self.syshape (useful for
        lattice graphs). If False, return as 1D array.
    backend : str, default 'scipy'
        Backend to use for computation: 'scipy', 'numpy', or 'cupy'.
    typf : type, default np.float64
        Data type for computation (e.g., np.float32, np.float64).
    transpose : bool, default True
        If True, stored as row-major (M×N) for indexing as self.eigV[i].
        If False, stored as column-major (N×M).
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    
    Returns
    -------
    ndarray
        The requested eigenvector with values in {-1, +1}, either as a 1D
        array or reshaped to system dimensions if reshaped=True.
    
    Notes
    -----
    This function is particularly useful when you need binary eigenvector
    representations for:
    - Community detection (sign indicates community membership)
    - State initialization in spin models
    - Feature extraction for classification
    - Visualization of eigenvector structure
    
    The binarization uses the `bin_sign()` function which maps:
    - Negative values → -1
    - Positive values → +1
    - Zero values → +1 (treating zero as positive)
    
    See Also
    --------
    get_eigV_check : More flexible version with binarization option
    get_eigV_binarized : Binarize without automatic computation
    bin_sign : The underlying binarization function
    
    Examples
    --------
    >>> # Get binarized eigenvector with automatic computation
    >>> v1_bin = get_eigV_bin_check(graph, which=1)
    >>> set(v1_bin)
    {-1, 1}
    
    >>> # Get binarized and reshaped eigenvector
    >>> v0_bin_grid = get_eigV_bin_check(graph, which=0, reshaped=True)
    >>> v0_bin_grid.shape
    (10, 10)  # For a 10x10 lattice
    >>> set(v0_bin_grid.ravel())
    {-1, 1}
    
    >>> # Use GPU acceleration
    >>> v2_bin = get_eigV_bin_check(graph, which=2, backend='cupy')
    """
    return get_eigV_check(
        self, which=which, binarize=True, reshaped=reshaped,
        backend=backend, typf=typf, transpose=transpose, flip_to_pos=flip_to_pos
    )

def get_eigV_check_list(
    self: "SignedGraph",
    custom_slice: slice = slice(None),
    binarize: bool = False,
    reshaped: bool = False,
    asarray: bool = False,
    backend: str = 'scipy',
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
) -> Union[List[NDArray], NDArray]:
    """
    Get multiple eigenvectors as a list or array with automatic computation.
    
    This function retrieves a range of eigenvectors specified by a slice object.
    It automatically computes eigenvectors if needed and returns them either as
    a list of arrays or as a single stacked array. Supports optional 
    binarization and reshaping.
    
    Parameters
    ----------
    custom_slice : slice, default slice(None)
        Slice object specifying which eigenvectors to retrieve. For example:
        - slice(None) retrieves all N eigenvectors (0 to N-1)
        - slice(5) retrieves eigenvectors 0 to 4
        - slice(2, 8) retrieves eigenvectors 2 to 7
        Note: Only the stop parameter is used; start and step are ignored.
    binarize : bool, default False
        If True, return binarized eigenvectors with values in {-1, +1}.
        If False, return raw continuous eigenvector values.
    reshaped : bool, default False
        If True, reshape eigenvectors to self.syshape (useful for lattice graphs).
        If False, return as 1D arrays (or 2D if asarray=True).
    asarray : bool, default False
        If True, return eigenvectors stacked as a 2D array with shape (k, N)
        where k is the number of eigenvectors and N is the number of nodes
        (or shape (k, *syshape) if reshaped=True).
        If False, return as a list of arrays.
    backend : str, default 'scipy'
        Backend to use for computation: 'scipy', 'numpy', or 'cupy'.
    typf : type, default np.float64
        Data type for computation (e.g., np.float32, np.float64).
    transpose : bool, default True
        If True, eigenvectors stored as row-major (M×N) for indexing as self.eigV[i].
        If False, stored as column-major (N×M).
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    
    Returns
    -------
    list of ndarray or ndarray
        If asarray=False: List of eigenvectors, each as an array.
        If asarray=True: Stacked array where each row (or first dimension)
        is an eigenvector.
        Values are continuous (binarize=False) or in {-1, +1} (binarize=True).
    
    Raises
    ------
    AssertionError
        If the requested slice stop index exceeds the number of nodes N.
    
    Notes
    -----
    This function optimizes computation by delegating to `get_eigV_check()`,
    which intelligently chooses between sparse (for small k) and dense
    (for large k) computation methods.
    
    The function always starts from eigenvector index 0 regardless of the
    slice's start parameter. The slice's step parameter is also ignored.
    
    See Also
    --------
    get_eigV_check : Get a single eigenvector with automatic computation
    get_eigV_bin_check : Get a single binarized eigenvector
    get_eigV_bin_check_list : Convenience wrapper for binarized eigenvectors
    get_sgspect_basis : Alternative for getting eigenvector basis
    compute_k_eigvV : Compute k smallest eigenvectors
    compute_laplacian_spectrum_weigV : Compute full spectrum
    
    Examples
    --------
    >>> # Get first 5 raw eigenvectors as a list
    >>> eigvecs = get_eigV_check_list(graph, custom_slice=slice(5))
    >>> len(eigvecs)
    5
    >>> eigvecs[0].shape
    (100,)
    
    >>> # Get first 5 binarized eigenvectors as a list
    >>> eigvecs_bin = get_eigV_check_list(
    ...     graph, custom_slice=slice(5), binarize=True
    ... )
    >>> set(eigvecs_bin[0])
    {-1, 1}
    
    >>> # Get all eigenvectors as a 2D array
    >>> eigvecs_arr = get_eigV_check_list(graph, asarray=True)
    >>> eigvecs_arr.shape
    (100, 100)  # For a graph with 100 nodes
    
    >>> # Get first 10 binarized eigenvectors as array, reshaped for lattice
    >>> eigvecs_10 = get_eigV_check_list(
    ...     graph, custom_slice=slice(10), binarize=True, 
    ...     reshaped=True, asarray=True
    ... )
    >>> eigvecs_10.shape
    (10, 10, 10)  # For a 10x10 lattice with 10 eigenvectors
    
    >>> # Use GPU acceleration for computation
    >>> eigvecs_gpu = get_eigV_check_list(
    ...     graph, custom_slice=slice(20), backend='cupy', asarray=True
    ... )
    """
    # Determine the maximum index to retrieve
    maxidx = custom_slice.stop if custom_slice.stop else self.N
    assert maxidx <= self.N, SG_ERRMSG_MAXEIGVIDX
    
    # Retrieve eigenvectors using the established method
    # Note: get_eigV_check already handles computation checks automatically
    eigvec_list = [
        get_eigV_check(
            self, i, binarize=binarize, reshaped=reshaped, backend=backend, 
            typf=typf, transpose=transpose, flip_to_pos=flip_to_pos
        ) 
        for i in range(maxidx)
    ]
    
    if asarray:
        return np.array(eigvec_list)
    return eigvec_list

def get_eigV_bin_check_list(
    self: "SignedGraph",
    custom_slice: slice = slice(None),
    reshaped: bool = False,
    asarray: bool = False,
    backend: str = 'scipy',
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
) -> Union[List[NDArray], NDArray]:
    """
    Get multiple binarized eigenvectors as a list or array.
    
    This is a convenience wrapper around `get_eigV_check_list()` that always
    returns binarized eigenvectors with values in {-1, +1}. It combines
    automatic computation with binarization for retrieving multiple eigenvectors.
    
    Parameters
    ----------
    custom_slice : slice, default slice(None)
        Slice object specifying which eigenvectors to retrieve. For example:
        - slice(None) retrieves all N eigenvectors (0 to N-1)
        - slice(5) retrieves eigenvectors 0 to 4
        - slice(2, 8) retrieves eigenvectors 2 to 7
        Note: Only the stop parameter is used; start and step are ignored.
    reshaped : bool, default False
        If True, reshape eigenvectors to self.syshape (useful for lattice graphs).
        If False, return as 1D arrays.
    asarray : bool, default False
        If True, return eigenvectors stacked as a 2D array with shape (k, N)
        where k is the number of eigenvectors and N is the number of nodes.
        If False, return as a list of 1D arrays.
    backend : str, default 'scipy'
        Backend to use for computation: 'scipy', 'numpy', or 'cupy'.
    typf : type, default np.float64
        Data type for computation (e.g., np.float32, np.float64).
    transpose : bool, default True
        If True, eigenvectors stored as row-major (M×N) for indexing as self.eigV[i].
        If False, stored as column-major (N×M).
    flip_to_pos : bool, default True
        If True, flip each eigenvector to have positive majority.
    
    Returns
    -------
    list of ndarray or ndarray
        If asarray=False: List of binarized eigenvectors, each with values
        in {-1, +1}.
        If asarray=True: 2D array where each row is a binarized eigenvector.
    
    Raises
    ------
    AssertionError
        If the requested slice stop index exceeds the number of nodes N.
    
    Notes
    -----
    This function is particularly useful when you need binary eigenvector
    representations for:
    - Community detection (sign indicates community membership)
    - State initialization in spin models
    - Feature extraction for classification
    - Batch processing of eigenvector structures
    
    See Also
    --------
    get_eigV_check_list : More flexible version with binarization option
    get_eigV_bin_check : Get a single binarized eigenvector
    
    Examples
    --------
    >>> # Get first 5 binarized eigenvectors as a list
    >>> eigvecs = get_eigV_bin_check_list(graph, custom_slice=slice(5))
    >>> len(eigvecs)
    5
    >>> set(eigvecs[0])
    {-1, 1}
    
    >>> # Get all binarized eigenvectors as an array
    >>> eigvecs_arr = get_eigV_bin_check_list(graph, asarray=True)
    >>> eigvecs_arr.shape
    (100, 100)  # For a graph with 100 nodes
    >>> np.all(np.isin(eigvecs_arr, [-1, 1]))
    True
    
    >>> # Use GPU acceleration
    >>> eigvecs_gpu = get_eigV_bin_check_list(
    ...     graph, custom_slice=slice(10), backend='cupy', asarray=True
    ... )
    """
    return get_eigV_check_list(
        self, custom_slice=custom_slice, binarize=True, 
        reshaped=reshaped, asarray=asarray,
        backend=backend, typf=typf, transpose=transpose, flip_to_pos=flip_to_pos
    )

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

#
# Spectral computations
# 






    

    

#
# Layout flip helpers
#
def make_eigV_transposed(self: "SignedGraph"):
    self.eigV = self.eigV.T
    self._eigV_is_transposed = True

def make_eigV_column_major(self: "SignedGraph"):
    self.eigV = self.eigV.T
    self._eigV_is_transposed = False


def make_adj_eigV_transposed(self: "SignedGraph"):
    self.adj_eigV = self.adj_eigV.T
    self._adj_eigV_is_transposed = True


def make_adj_eigV_column_major(self: "SignedGraph"):
    self.adj_eigV = self.adj_eigV.T
    self._adj_eigV_is_transposed = False

    
