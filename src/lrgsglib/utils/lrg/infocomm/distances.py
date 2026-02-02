"""Distance functions based on Laplacian spectrum.

This module provides functions to compute spectral distances between graph nodes.
"""

import numpy as np
from numpy.typing import NDArray
from scipy.cluster.hierarchy import cophenet
from scipy.linalg import expm
from scipy.sparse import csr_matrix, csr_array
from scipy.sparse.linalg import expm as sparse_expm
from scipy.spatial.distance import squareform

__all__ = [
    "extract_ultrametric_matrix",
    "lapl_dists",
]


def extract_ultrametric_matrix(linkage_matrix: NDArray, n_nodes: int) -> NDArray:
    """
    Extract the ultrametric distance matrix from a linkage matrix.

    Parameters
    ----------
    linkage_matrix : np.ndarray
        The linkage matrix from hierarchical clustering.
    n_nodes : int
        Number of original nodes/observations.

    Returns
    -------
    np.ndarray
        The ultrametric distance matrix (symmetric, n_nodes x n_nodes).
    """
    # Compute cophenetic distances (these are the ultrametric distances)
    cophenetic_dists = cophenet(linkage_matrix)

    # Convert from condensed form to square matrix
    ultrametric_matrix = squareform(cophenetic_dists)

    return ultrametric_matrix


def lapl_dists(L, tau: float = 1e-2, is_signed: bool = False) -> NDArray:
    """
    Compute the pairwise distances between nodes based on the Laplacian spectrum.

    Parameters
    ----------
    L : NDArray or csr_matrix
        The Laplacian matrix of the graph. Can be a dense numpy array or a
        sparse scipy CSR matrix.
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
            dists = np.sqrt(1 - old_d / np.max(old_d))
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
