import numpy as np
#
from typing import Tuple, Callable
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
    "entropy"
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

    # Generate logarithmically spaced time steps
    t = np.logspace(t1, t2, steps)

    # Initialize arrays for entropy and variance
    S = np.zeros(len(t))
    VarL = np.zeros(len(t))

    # Compute entropy and variance for each time step
    for i, tau in enumerate(t):
        rhoTr = np.exp(-tau * w)  # Exponential decay of eigenvalues
        Tr = np.nansum(rhoTr)    # Trace of the exponential matrix
        rho = rhoTr / Tr         # Normalized eigenvalue distribution

        # Compute normalized entropy
        S[i] = -np.nansum(rho * np.log(rho)) / np.log(N)

        # Compute variance of the Laplacian spectrum
        avgrho = np.nansum(w * rhoTr) / Tr
        av2rho = np.nansum((w**2) * rhoTr) / Tr
        VarL[i] = av2rho - avgrho**2

    # Compute the derivative of entropy with respect to time
    dS = np.log(N) * np.diff(1 - S) / np.diff(np.log(t))

    return 1 - S, dS, VarL, t