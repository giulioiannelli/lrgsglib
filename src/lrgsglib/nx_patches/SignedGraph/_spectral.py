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
def compute_laplacian_spectrum(
    self: "SignedGraph",
    typf: type = np.float64,
    backend: Optional[str] = None,
    keep_sparse: bool = None,
    verbose: bool = False,
):
    """
    Compute eigenvalues only of the signed Laplacian (no eigenvectors).

    This function computes only the eigenvalues of the signed Laplacian matrix
    without computing eigenvectors, which is faster and uses less memory when
    eigenvectors are not needed.

    Parameters
    ----------
    typf : type, default np.float64
        Data type for computation (e.g., np.float32 for memory savings,
        np.float64 for higher precision).
    backend : str, optional
        Override the instance backend ('numpy', 'scipy', or 'cupy').
        Defaults to the backend selected on the graph instance.
    keep_sparse : bool, optional
        If True, use sparse methods (computes N-2 eigenvalues).
        If False, convert to dense.
        If None (default), automatically decide based on N and sparsity.

    Returns
    -------
    None
        Results are stored in instance attribute:

        - self.eigv : ndarray of shape (N,) or (N-2,)
            Eigenvalues in ascending order. If sparse method used, returns N-2.

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

    **Performance:**
    - Dense: Uses `eigvalsh` (eigenvalues-only, faster than `eigh`)
    - Sparse: Uses `eigsh` (gets N-2 eigenvalues, discards eigenvectors)

    **Use Cases:**
    - Computing spectral invariants (trace, determinant, spectral radius)
    - Analyzing spectral gap
    - Studying energy landscapes
    - Entropy calculations (N-2 is sufficient for large N)
    - Memory-constrained environments

    **Sparse Eigenvalue Limitation:**
    Sparse eigensolvers (ARPACK) typically return N-2 eigenvalues due to
    numerical limitations. This is stored in the ``_sparse_spectrum`` flag.
    For full N eigenvalues, use dense computation (``keep_sparse=False``) or
    ensure the graph is small enough for dense methods. The eigenvalues are
    stored in ``self.eigv``. Use ``self._sparse_spectrum`` to check if the
    spectrum is complete (False) or sparse (True, N-2 values).

    See Also
    --------
    compute_laplacian_spectrum_weigV : Compute both eigenvalues and eigenvectors
    compute_k_eigvV : Compute k smallest eigenvalues and eigenvectors

    Examples
    --------
    >>> # Compute eigenvalues only (auto-selects dense/sparse)
    >>> compute_laplacian_spectrum(graph)
    >>> graph.eigv.shape
    (100,)

    >>> # Force sparse method for large graph
    >>> compute_laplacian_spectrum(graph, keep_sparse=True)
    >>> graph.eigv.shape
    (9998,)  # N-2 for N=10000

    >>> # Access spectral properties
    >>> spectral_gap = graph.eigv[1] - graph.eigv[0]
    >>> trace = graph.eigv.sum()
    """
    from ._backend import BackendManager, Backend

    # Check if we have a valid cached spectrum
    cached_eigv = getattr(self, "eigv", None)
    has_cache = cached_eigv is not None and len(cached_eigv) >= self.N - 2

    if has_cache and keep_sparse is None:
        # Auto-mode: accept cached spectrum (sparse or dense)
        return

    if has_cache and keep_sparse is not None:
        # Explicit mode: check if cached spectrum matches requested sparsity
        cached_is_sparse = len(cached_eigv) == self.N - 2
        cached_is_dense = len(cached_eigv) == self.N

        if keep_sparse and cached_is_sparse:
            # Requested sparse, have sparse -> use cache
            return
        elif not keep_sparse and cached_is_dense:
            # Requested dense, have dense -> use cache
            return
        # Otherwise fall through to recompute

    backend_requested = backend
    if backend is None:
        # Use instance default if available, fall back to numpy
        backend = getattr(self, "_backend_name", Backend.NUMPY.value)

    # Always go through BackendManager for consistent validation and fallback
    backend_obj = BackendManager.get_backend(backend, fallback=True)

    # Ensure signed Laplacian matrix is computed
    if self.slp is None:
        self.upd_graph_matrices()

    # Determine if we should use sparse methods
    if keep_sparse is None:
        # Auto-decide based on size and sparsity
        import scipy.sparse

        # Detect sparsity more robustly
        if hasattr(self.slp, 'nnz'):
            sparsity = self.slp.nnz / (self.N * self.N)
        elif scipy.sparse.issparse(self.slp):
            sparsity = self.slp.count_nonzero() / (self.N * self.N)
        else:
            # Dense matrix
            sparsity = 1.0

        # For CuPy backend, always use dense (full spectrum needed for entropy)
        if backend_obj.name == 'cupy':
            keep_sparse = False
        else:
            keep_sparse = self.N > 5000 and sparsity < 0.3

    self._sparse_spectrum = keep_sparse
    self._last_spectrum_backend = backend_obj.name
    self._last_spectrum_request = backend_requested
    self.spectrum_params = {
        "backend_requested": backend_requested,
        "backend_used": backend_obj.name,
        "keep_sparse": keep_sparse,
        "with_eigvectors": False,
    }

    # Memory tracking for debugging
    if verbose:
        try:
            import psutil
            import os
            process = psutil.Process(os.getpid())
            mem_before_spectrum = process.memory_info().rss / (1024**3)
            print(f"[MEM-SPECTRAL] Before spectrum computation: {mem_before_spectrum:.2f} GB", flush=True)
            print(f"[INFO-SPECTRAL] N={self.N:,}, backend={backend_obj.name}, keep_sparse={keep_sparse}", flush=True)
        except:
            pass

    if keep_sparse:
        # Use sparse methods - gets N-2 eigenvalues, discards eigenvectors
        # Cast to correct dtype only if needed
        laplacian = self.slp if self.slp.dtype == typf else self.slp.astype(typf)
        self.eigv, _ = backend_obj.eigh_sparse(
            laplacian,
            k=None,
            return_eigenvectors=False,
        )
    else:
        # Convert to dense for full eigenvalue computation
        # CRITICAL: Avoid intermediate copies to minimize peak memory
        if verbose:
            try:
                print(f"[STEP 3/5] Converting sparse Laplacian to dense ({self.N}×{self.N})...", flush=True)
                sparse_mb = (self.slp.nnz * 16) / (1024**2) if hasattr(self.slp, 'nnz') else 0
                dense_mb = (self.N ** 2 * 8) / (1024**2)
                print(f"[INFO-SPECTRAL] Sparse: {sparse_mb:.1f} MB → Dense: {dense_mb:.1f} MB", flush=True)
            except:
                pass

        if hasattr(self.slp, "toarray"):
            # Sparse matrix: check dtype before converting
            if self.slp.dtype == typf:
                # Already correct dtype: direct conversion
                dense_laplacian = self.slp.toarray()
            else:
                # Wrong dtype: must copy during conversion
                # Cast sparse first (cheap), then convert to dense
                dense_laplacian = self.slp.astype(typf).toarray()

            if verbose:
                try:
                    mem_after_dense = process.memory_info().rss / (1024**3)
                    print(f"[MEM-SPECTRAL] After dense conversion: {mem_after_dense:.2f} GB (+{mem_after_dense - mem_before_spectrum:.2f} GB)", flush=True)
                except:
                    pass

            # CRITICAL: Explicitly clear sparse Laplacian to free memory before eigval computation
            # For very large matrices (N>100k), this can save 10-50GB of heap space
            # The sparse matrix is no longer needed after dense conversion
            sparse_key = self.on_g
            if sparse_key in self.signed_laplacian_matrices:
                del self.signed_laplacian_matrices[sparse_key]
                if verbose:
                    print(f"[INFO-SPECTRAL] Deleted sparse Laplacian from cache", flush=True)
            import gc
            gc.collect()  # Force garbage collection to reclaim sparse matrix memory

            if verbose:
                try:
                    mem_after_gc = process.memory_info().rss / (1024**3)
                    print(f"[MEM-SPECTRAL] After GC: {mem_after_gc:.2f} GB (freed {mem_after_dense - mem_after_gc:.2f} GB)", flush=True)
                except:
                    pass
        else:
            # Already dense: just cast if needed
            dense_laplacian = np.asarray(self.slp, dtype=typf)

        # Use eigvalsh (eigenvalues only, faster than eigh)
        if verbose:
            print(f"[STEP 4/5] Computing eigenvalues with {backend_obj.name}.eigvalsh...", flush=True)

        self.eigv = backend_obj.eigvalsh(dense_laplacian)

        if verbose:
            try:
                mem_after_eigvalsh = process.memory_info().rss / (1024**3)
                print(f"[MEM-SPECTRAL] After eigvalsh: {mem_after_eigvalsh:.2f} GB", flush=True)
                print(f"[STEP 5/5] Eigenvalue computation complete: {len(self.eigv)} values", flush=True)
            except:
                pass
#
def compute_k_eigvV(
    self: "SignedGraph",
    k: int = 1,
    backend: Optional[str] = None,
    transpose: bool = True,
    flip_to_pos: bool = True,
    typf: type = np.float64,
    which: str = 'SM',
    solver: str = 'eigsh',
):
    """
    Compute k smallest eigenvectors of the signed Laplacian.
    
    Parameters
    ----------
    k : int, default 1
        Number of eigenvectors to compute.
    backend : str, default 'scipy'
        Backend to use: 'scipy', 'numpy', or 'cupy'.
    solver : str, default 'eigsh'
        Sparse solver to use when applicable: 'eigsh' (ARPACK) or 'lobpcg'.
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
    backend_name = backend if backend is not None else getattr(self, "_backend_name", "scipy")
    solver = solver.lower()
    if solver not in {"eigsh", "lobpcg"}:
        raise ValueError("solver must be 'eigsh' or 'lobpcg'.")

    if solver == "lobpcg":
        import warnings
        if k > self.N // 2:
            raise ValueError("lobpcg is intended for small k (k <= N/2).")

        if which in {"SM", "SA"}:
            if which == "SM":
                warnings.warn(
                    "lobpcg uses algebraic ordering; treating which='SM' as 'SA'.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            largest = False
        elif which in {"LM", "LA"}:
            if which == "LM":
                warnings.warn(
                    "lobpcg uses algebraic ordering; treating which='LM' as 'LA'.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            largest = True
        else:
            raise ValueError("lobpcg supports which in {'SM','SA','LM','LA'}.")

        if backend_name == "cupy":
            import cupy as cp
            import cupyx.scipy.sparse as cusp
            import cupyx.scipy.sparse.linalg as cusp_linalg
            from scipy.sparse import issparse, csr_matrix

            laplacian = self.slp.astype(typf)
            if issparse(laplacian):
                laplacian_gpu = cusp.csr_matrix(laplacian)
            else:
                laplacian_gpu = cusp.csr_matrix(csr_matrix(laplacian))

            X = cp.eye(self.N, k, dtype=typf)
            eigvals_gpu, eigvecs_gpu = cusp_linalg.lobpcg(
                laplacian_gpu, X, largest=largest
            )
            eigvals = cp.asnumpy(eigvals_gpu)
            eigvecs = cp.asnumpy(eigvecs_gpu)
        elif backend_name.startswith("scipy"):
            import scipy.sparse as scsp
            from scipy.sparse.linalg import lobpcg

            laplacian = self.slp.astype(typf)
            if not scsp.issparse(laplacian):
                laplacian = scsp.csr_matrix(laplacian)

            X = np.eye(self.N, k, dtype=typf)
            eigvals, eigvecs = lobpcg(laplacian, X, largest=largest)
        else:
            raise ValueError("lobpcg is supported only for scipy or cupy backends.")

        order = np.argsort(eigvals)
        if largest:
            order = order[::-1]
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]

        if flip_to_pos:
            eigvecs = flip_to_positive_majority_adapted(eigvecs, axis=0)

        self.eigv, self.eigV = eigvals, eigvecs
        self._eigV_is_transposed = False
        if transpose:
            make_eigV_transposed(self)
        return

    if (backend_name in ['numpy', 'cupy']) or k > self.N // 2:
        compute_laplacian_spectrum_weigV(
            self, backend_name, transpose, flip_to_pos, typf
        )
    elif backend_name.startswith('scipy'):
        mode = backend_name.split('_')
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
    backend: Optional[str] = None,
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
    backend_name = backend if backend is not None else getattr(self, "_backend_name", "scipy")
    if (backend_name in ['numpy', 'cupy']) or k > self.N // 2:
        compute_adjacency_spectrum_weigV(
            self, backend=backend_name, transpose=transpose,
            flip_to_pos=flip_to_pos, typf=typf
        )
        return

    if not backend_name.startswith('scipy'):
        raise ValueError(
            f"Unsupported backend '{backend_name}' for sparse adjacency eigensolver. "
            "Use a scipy backend or fall back to dense computation."
        )

    mode = backend_name.split('_')
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
    backend: Optional[str] = None,
    transpose: bool = True,
    flip_to_pos: bool = True,
    typf: type = np.float64,
    keep_sparse: bool = None,
):
    """
    Compute full eigendecomposition of the signed Laplacian.
    
    This function computes all eigenvalues and eigenvectors of the signed
    Laplacian matrix. It automatically chooses between dense and sparse methods
    based on matrix size and sparsity for optimal performance.
    
    Parameters
    ----------
    backend : str, optional
        Backend to use for computation. If None (default), uses self._backend
        set during initialization. Supported: 'numpy', 'scipy', 'cupy'.
        For backward compatibility, you can override the instance backend.
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
    keep_sparse : bool, optional
        If True, use sparse eigensolvers even for full spectrum (computes N-2 eigenpairs).
        If False, convert to dense and use dense methods.
        If None (default), automatically decide based on N and sparsity:
        - Use sparse for N > 5000 and sparsity < 0.3
        - Use dense otherwise
    
    Returns
    -------
    None
        Results are stored in instance attributes:

        - self.eigv : ndarray of shape (N,) or (N-2,) if sparse
            Eigenvalues in ascending order.
        - self.eigV : ndarray of shape (M, N) or (N, M)
            Eigenvectors corresponding to eigenvalues. Shape depends on
            transpose parameter.
        - self._eigV_is_transposed : bool
            Flag indicating current storage layout.
        - self._sparse_spectrum : bool
            Flag indicating if sparse method was used (N-2 eigenpairs)
    
    Notes
    -----
    **Performance Improvements:**
    
    - For N < 5000: Uses dense methods with divide-and-conquer (~2x faster)
    - For N > 5000 with sparse matrices: Uses iterative sparse solvers
      (saves memory, O(N*k) vs O(N²))
    
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
    - Use 'scipy' for best performance on CPU (uses optimized LAPACK drivers)
    - Use 'cupy' for large graphs where GPU acceleration provides significant
      speedup (requires cupy installation and compatible GPU)
    
    See Also
    --------
    compute_k_eigvV : Compute k smallest eigenvectors (efficient for small k)
    make_eigV_transposed : Convert to row-major storage
    make_eigV_column_major : Convert to column-major storage
    
    Examples
    --------
    >>> # Compute full spectrum with default settings (auto-selects method)
    >>> compute_laplacian_spectrum_weigV(graph)
    >>> graph.eigv.shape
    (100,)
    >>> graph.eigV.shape
    (100, 100)
    
    >>> # Force sparse method for large systems
    >>> compute_laplacian_spectrum_weigV(graph, keep_sparse=True)
    >>> graph.eigv.shape
    (9998,)  # N-2 eigenvalues for N=10000
    
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
        and len(cached_eigv) >= self.N - 2  # Accept N or N-2
        and cached_eigV.shape[0] >= self.N - 2
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

    # Use instance backend if not explicitly specified
    backend_requested = backend
    if backend is None:
        backend_obj = self._backend
    else:
        # Allow override for backward compatibility
        from ._backend import BackendManager
        backend_obj = BackendManager.get_backend(backend, fallback=True)

    # Ensure signed Laplacian matrix is computed
    if self.slp is None:
        self.upd_graph_matrices()

    # Determine if we should use sparse methods
    if keep_sparse is None:
        # Auto-decide based on size and sparsity
        sparsity = self.slp.nnz / (self.N * self.N) if hasattr(self.slp, 'nnz') else 1.0
        # For CuPy backend, always use dense (full spectrum needed for entropy)
        if backend_obj.name == 'cupy':
            keep_sparse = False
        else:
            keep_sparse = self.N > 5000 and sparsity < 0.3
    
    self._sparse_spectrum = keep_sparse
    self._last_spectrum_backend = backend_obj.name
    self._last_spectrum_request = backend_requested
    self.spectrum_params = {
        "backend_requested": backend_requested,
        "backend_used": backend_obj.name,
        "keep_sparse": keep_sparse,
        "with_eigvectors": True,
    }
    
    if keep_sparse:
        # Use sparse iterative methods (computes N-2 eigenpairs)
        logger.info(f"Using sparse methods for large system (N={self.N})")
        if backend_obj.name == 'cupy' and self.N < 10000:
            logger.warning(
                f"CuPy with sparse mode for N={self.N}. Consider using dense "
                "(keep_sparse=False) or scipy backend for better performance."
            )
        laplacian = self.slp.astype(typf)
        self.eigv, self.eigV = backend_obj.eigh_sparse(laplacian, k=None)
    else:
        # Convert sparse matrix to dense array for eigendecomposition
        slp = self.slp.astype(typf).toarray()
        
        # Compute eigenvalues and eigenvectors using backend
        self.eigv, self.eigV = backend_obj.eigh(slp)

    if flip_to_pos:
        self.eigV = flip_to_positive_majority_adapted(self.eigV, axis=0)

    self._eigV_is_transposed = False
    if transpose:
        make_eigV_transposed(self)


def compute_adjacency_spectrum_weigV(
    self: "SignedGraph",
    backend: Optional[str] = None,
    transpose: bool = True,
    flip_to_pos: bool = True,
    typf: type = np.float64
):
    """
    Compute the full eigendecomposition of the adjacency matrix.

    Parameters
    ----------
    backend : str, optional
        Backend to use for computation. If None (default), uses self._backend.
        Supported: 'numpy', 'scipy', 'cupy'.
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

    # Use instance backend if not explicitly specified
    if backend is None:
        backend_obj = self._backend
    else:
        # Allow override for backward compatibility
        from ._backend import BackendManager
        backend_obj = BackendManager.get_backend(backend, fallback=True)

    # Convert adjacency matrix to dense array
    adj_matrix = self.adj.astype(typf)
    dense_adj = adj_matrix.toarray() if hasattr(adj_matrix, "toarray") else np.asarray(adj_matrix, dtype=typf)

    # Compute eigenvalues and eigenvectors using backend
    self.adj_eigv, self.adj_eigV = backend_obj.eigh(dense_adj)

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
    backend: Optional[str] = None,
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
    backend_name = backend if backend is not None else getattr(self, "_backend_name", "scipy")
    if not hasattr(self, "eigV") or which >= len(self.eigV):
        if which < self.N // 2:
            compute_k_eigvV(
                self, k=which + 1, backend=backend_name, typf=typf,
                transpose=transpose, flip_to_pos=flip_to_pos
            )
        else:
            compute_laplacian_spectrum_weigV(
                self, typf=typf, transpose=transpose,
                backend=backend_name, flip_to_pos=flip_to_pos
            )
    #
    return get_eigV(self, which=which, binarize=binarize, reshape=reshaped)


def get_eigV_bin_check(
    self: "SignedGraph",
    which: int = 0,
    reshaped: bool = False,
    backend: Optional[str] = None,
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
    backend_name = backend if backend is not None else getattr(self, "_backend_name", "scipy")
    return get_eigV_check(
        self, which=which, binarize=True, reshaped=reshaped,
        backend=backend_name, typf=typf, transpose=transpose, flip_to_pos=flip_to_pos
    )

def get_eigV_check_list(
    self: "SignedGraph",
    custom_slice: slice = slice(None),
    binarize: bool = False,
    reshaped: bool = False,
    asarray: bool = False,
    backend: Optional[str] = None,
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
    backend_name = backend if backend is not None else getattr(self, "_backend_name", "scipy")
    eigvec_list = [
        get_eigV_check(
            self, i, binarize=binarize, reshaped=reshaped, backend=backend_name, 
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
    backend: Optional[str] = None,
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
    backend_name = backend if backend is not None else getattr(self, "_backend_name", "scipy")
    return get_eigV_check_list(
        self, custom_slice=custom_slice, binarize=True, 
        reshaped=reshaped, asarray=asarray,
        backend=backend_name, typf=typf, transpose=transpose, flip_to_pos=flip_to_pos
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


#
# Quantum propagator methods
#
def compute_quantum_propagator(
    self: "SignedGraph",
    t: float,
    full_matrix: bool = False
) -> Union[tuple[NDArray, NDArray], NDArray]:
    """
    Compute quantum propagator U(t) = exp(-i*t*L) for this graph.

    Parameters
    ----------
    t : float
        Quantum time parameter.
    full_matrix : bool, optional
        If True, return full U(t) matrix. If False, return phases and
        eigenvectors (default: False).

    Returns
    -------
    If full_matrix=False:
        phases : NDArray
            Quantum phases exp(-i*t*λ)
        eigenvectors : NDArray
            Laplacian eigenvectors
    If full_matrix=True:
        U_t : NDArray
            Full quantum propagator matrix

    Raises
    ------
    ValueError
        If eigendecomposition has not been computed.

    Notes
    -----
    Requires eigendecomposition. Call compute_k_eigvV() or
    compute_laplacian_spectrum_weigV() first.

    The quantum propagator is unitary: U(t)^† U(t) = I

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> l = Lattice2D(side=10, pflip=0.1)
    >>> l.flip_random_fract_edges()
    >>> l.compute_k_eigvV()
    >>> U_t = l.compute_quantum_propagator(t=1.0, full_matrix=True)
    >>> print(np.allclose(U_t @ U_t.conj().T, np.eye(l.N)))
    True
    """
    if not hasattr(self, 'eigv') or not hasattr(self, 'eigV'):
        raise ValueError(
            "Eigendecomposition required. "
            "Call compute_k_eigvV() or compute_laplacian_spectrum_weigV() first."
        )

    from ...utils.lrg.quantum import compute_quantum_propagator_spectral
    return compute_quantum_propagator_spectral(
        self.eigv, self.eigV, t, full_matrix=full_matrix
    )


def quantum_walk_probabilities(
    self: "SignedGraph",
    t: float,
    init_node: int
) -> NDArray:
    """
    Compute quantum walk probability distribution at time t.

    Parameters
    ----------
    t : float
        Quantum time.
    init_node : int
        Initial node (source of quantum walk).

    Returns
    -------
    probabilities : NDArray, shape (N,)
        Quantum probability at each node: P_j(t) = \\|⟨j\\|U(t)\\|i⟩\\|²

    Raises
    ------
    ValueError
        If eigendecomposition has not been computed or init_node is invalid.

    Notes
    -----
    Starting from localized initial state at init_node.
    Probability is conserved: sum(probabilities) = 1.0

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> l = Lattice2D(side=10)
    >>> l.compute_k_eigvV()
    >>> P_t = l.quantum_walk_probabilities(t=1.0, init_node=50)
    >>> print(np.isclose(np.sum(P_t), 1.0))
    True
    """
    if not hasattr(self, 'eigv') or not hasattr(self, 'eigV'):
        raise ValueError(
            "Eigendecomposition required. "
            "Call compute_k_eigvV() or compute_laplacian_spectrum_weigV() first."
        )

    if init_node < 0 or init_node >= self.N:
        raise ValueError(f"init_node={init_node} out of range [0, {self.N})")

    from ...utils.lrg.quantum import quantum_probability_distribution
    return quantum_probability_distribution(
        self.eigv, self.eigV, t, init_node
    )


def quantum_observables_time_series(
    self: "SignedGraph",
    t_array: NDArray,
    init_type: str = "localized",
    init_node: int = 0,
) -> dict[str, NDArray]:
    """
    Compute quantum observables over time array.

    Parameters
    ----------
    t_array : NDArray
        Array of quantum time values.
    init_type : str, optional
        Initial state type: "uniform", "localized", "thermal"
        (default: "localized").
    init_node : int, optional
        Initial node for localized state (default: 0).

    Returns
    -------
    observables : dict
        Dictionary with keys:
        - "time": Time array
        - "von_neumann_entropy": von Neumann entropy S_vN(t)
        - "coherence": Quantum coherence C(t)
        - "purity": Purity Tr(ρ²)
        - "participation_ratio": Inverse participation ratio
        - "probability_distribution": P_j(t) for all nodes

    Raises
    ------
    ValueError
        If eigendecomposition has not been computed.

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> import numpy as np
    >>> l = Lattice2D(side=10)
    >>> l.compute_k_eigvV()
    >>> t_vals = np.linspace(0, 10, 100)
    >>> obs = l.quantum_observables_time_series(t_vals, init_node=50)
    >>> print(obs.keys())
    dict_keys(['time', 'von_neumann_entropy', 'coherence', 'purity',
               'participation_ratio', 'probability_distribution'])
    """
    if not hasattr(self, 'eigv') or not hasattr(self, 'eigV'):
        raise ValueError(
            "Eigendecomposition required. "
            "Call compute_k_eigvV() or compute_laplacian_spectrum_weigV() first."
        )

    from ...utils.lrg.quantum import compute_quantum_observables_from_eigenvalues
    return compute_quantum_observables_from_eigenvalues(
        self.eigv, self.eigV, self.N, t_array, init_type, init_node
    )
