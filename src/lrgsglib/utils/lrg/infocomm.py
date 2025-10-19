import numpy as np
#
from typing import Tuple, Callable, Optional
from networkx import Graph
#
from numpy.typing import NDArray
from scipy.cluster.hierarchy import cophenet
from scipy.linalg import expm
from scipy.sparse import csr_matrix, csr_array
from scipy.sparse.linalg import expm as sparse_expm
from scipy.spatial.distance import squareform
#
from ..basic import dtype_numerical_precision
from .spectral import *
#
__all__ = [
    "extract_ultrametric_matrix",
    "lapl_dists",
    "entropy",
    "compute_entropy_observables_from_eigenvalues",
]
#
def extract_ultrametric_matrix(linkage_matrix, n_nodes):
    """
    Extract the ultrametric distance matrix from a linkage matrix.
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        The linkage matrix from hierarchical clustering
    n_nodes : int
        Number of original nodes/observations
        
    Returns:
    --------
    np.ndarray
        The ultrametric distance matrix (symmetric, n_nodes x n_nodes)
    """
    # Compute cophenetic distances (these are the ultrametric distances)
    cophenetic_dists = cophenet(linkage_matrix)
    
    # Convert from condensed form to square matrix
    ultrametric_matrix = squareform(cophenetic_dists)
    
    return ultrametric_matrix
#
def lapl_dists(L, tau: float = 1e-2, is_signed: bool = False) -> NDArray:
    """
    Compute the pairwise distances between nodes based on the Laplacian spectrum.

    Parameters
    ----------
    L : NDArray or csr_matrix
        The Laplacian matrix of the graph. Can be a dense numpy array or a sparse scipy CSR matrix.
    tau : float, optional
        A scaling parameter for the exponential function (default is 1e-2).
    is_signed : bool, optional
        Whether the graph is signed (default is False).

    Returns
    -------
    NDArray
        A condensed distance matrix representing pairwise distances between nodes.
    """
    if isinstance(L, (csr_matrix, csr_array)):
        # Sparse matrix computation
        num = sparse_expm(-tau * L)
        trace_num = num.diagonal().sum()
        rho = num / trace_num
        Trho = rho.copy()
        Trho.data = 1.0 / rho.data
        Trho = Trho.maximum(Trho.T)
        Trho.setdiag(0)
        if is_signed:
            old_d = squareform(Trho.toarray())
            dists = np.sqrt(1 - old_d/np.max(old_d))
        else:
            dists = squareform(Trho.toarray())
    else:
        # Dense matrix computation
        num = expm((-tau * L))
        rho = num / np.trace(num)
        Trho = np.copy(1.0 / rho)
        Trho = np.maximum(Trho, Trho.T)
        np.fill_diagonal(Trho, 0)
        if is_signed:
            old_d = squareform(Trho)
            dists = np.sqrt(np.max(old_d) - old_d)
        else:
            dists = squareform(Trho)
    return dists


def entropy(
    G: Graph,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    wTresh: float = dtype_numerical_precision(np.float64),
    func_lspectrum: Callable = get_graph_lspectrum
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Compute the entropy, its derivative, variance, and time steps for a graph.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    steps : int, optional
        Number of time steps for the computation (default is 600).
    t1 : int, optional
        Logarithmic start time (default is -2).
    t2 : int, optional
        Logarithmic end time (default is 5).
    wTresh : float, optional
        Threshold for filtering eigenvalues (default is machine precision for float64).
    func_lspectrum : Callable, optional
        Function to compute the Laplacian spectrum of the graph (default is `get_graph_lspectrum`).

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray, NDArray]
        - Normalized entropy values.
        - Derivative of entropy with respect to time.
        - Variance of the Laplacian spectrum.
        - Time steps used for the computation.
    """
    # Number of nodes in the graph
    N = G.number_of_nodes()

    # Compute the Laplacian spectrum
    L, w = func_lspectrum(G)

    # Filter eigenvalues based on the threshold
    w = np.where(np.abs(w) > wTresh, w, 0)

    normalized_entropy, entropy_derivative, variance_profile, time_grid = (
        compute_entropy_observables_from_eigenvalues(
            eigenvalues=w,
            num_nodes=N,
            steps=steps,
            t1=t1,
            t2=t2,
            typf=np.float64,
            threshold=wTresh,
            pad_last=False,
        )
    )

    return normalized_entropy, entropy_derivative, variance_profile, time_grid


def compute_entropy_observables_from_eigenvalues(
    eigenvalues: NDArray,
    num_nodes: int,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    typf: type = np.float64,
    threshold: Optional[float] = None,
    pad_last: bool = False,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Compute entropy-related observables from a Laplacian spectrum.

    Parameters
    ----------
    eigenvalues : NDArray
        Spectrum of the Laplacian (or signed Laplacian) matrix.
    num_nodes : int
        Number of nodes in the graph used to normalize the entropy.
    steps : int, optional
        Number of logarithmic time samples (default 600).
    t1 : int, optional
        Lower exponent for the log-spaced time grid (default -2).
    t2 : int, optional
        Upper exponent for the log-spaced time grid (default 5).
    typf : type, optional
        Floating-point type used for computations (default np.float64).
    threshold : float, optional
        Magnitude threshold below which eigenvalues are treated as zero.
        Defaults to the machine precision of `typf`.
    pad_last : bool, optional
        When True, append a final derivative sample equal to the last computed
        value so that the derivative array matches the time grid length.

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray, NDArray]
        Normalized entropy profile (1 - S), entropy derivative, variance profile,
        and the sampled time grid.
    """
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    eigvals = np.asarray(eigenvalues, dtype=typf)
    eps = threshold if threshold is not None else dtype_numerical_precision(typf)
    eigvals = np.where(np.abs(eigvals) > eps, eigvals, typf(0))

    time_grid = np.logspace(t1, t2, steps, dtype=typf)
    entropy_profile = np.zeros_like(time_grid, dtype=typf)
    variance_profile = np.zeros_like(time_grid, dtype=typf)

    log_N = np.log(typf(num_nodes)) if num_nodes > 0 else typf(1)

    for idx, tau in enumerate(time_grid):
        rhoTr = np.exp(-tau * eigvals)
        trace_rho = np.nansum(rhoTr, dtype=typf)
        if trace_rho == 0:
            rho = np.zeros_like(rhoTr)
        else:
            rho = rhoTr / trace_rho

        with np.errstate(divide='ignore', invalid='ignore'):
            entropy_profile[idx] = -np.nansum(rho * np.log(rho), dtype=typf) / log_N

        if trace_rho:
            avgrho = np.nansum(eigvals * rhoTr, dtype=typf) / trace_rho
            av2rho = np.nansum((eigvals ** 2) * rhoTr, dtype=typf) / trace_rho
        else:
            avgrho = typf(0)
            av2rho = typf(0)
        variance_profile[idx] = av2rho - avgrho ** 2

    normalized_entropy = 1 - entropy_profile

    if steps > 1:
        entropy_derivative_core = np.asarray(
            log_N * np.diff(normalized_entropy) / np.diff(np.log(time_grid)),
            dtype=typf,
        )
        if pad_last:
            entropy_derivative = np.empty_like(time_grid, dtype=typf)
            entropy_derivative[:-1] = entropy_derivative_core
            entropy_derivative[-1] = entropy_derivative_core[-1]
        else:
            entropy_derivative = entropy_derivative_core
    else:
        entropy_derivative = (
            np.zeros_like(time_grid, dtype=typf) if pad_last else np.empty(0, dtype=typf)
        )

    return normalized_entropy, entropy_derivative, variance_profile, time_grid
