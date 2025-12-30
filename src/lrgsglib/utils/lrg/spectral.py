import networkx as nx
import numpy as np
#
from typing import Optional
#
from numpy.typing import NDArray
from scipy.linalg import fractional_matrix_power
from scipy.linalg import expm
from scipy.sparse import diags
#
__all__ = [
    "get_graph_lspectrum",
    "get_graph_lspectrum_rw",
    "compute_laplacian_properties",
]

def get_graph_lspectrum(
    G: nx.Graph,
    library: str = "numpy",
) -> tuple[NDArray, NDArray]:
    """
    Compute the Laplacian matrix and its spectrum for a given graph.

    Parameters
    ----------
    G : nx.Graph
        The input graph for which the Laplacian matrix and spectrum are computed.
    library : str, optional
        The library to use for eigenvalue computation. Options include
        ``"numpy"``, ``"networkx"``/``"nx"``, ``"scipy"``/``"sp"``,
        and ``"cupy"``/``"cp"``. Defaults to ``"numpy"``.

    Returns
    -------
    tuple[NDArray, NDArray]
        A tuple containing:

        - L : NDArray
            The Laplacian matrix of the graph.
        - w : NDArray
            The spectrum (eigenvalues) of the Laplacian matrix. Eigenvalues
            are not sorted.

    Notes
    -----
    - The Laplacian matrix is computed using NetworkX's `laplacian_matrix`.
    - The spectrum is computed using the specified library for eigenvalue computation.
    - If ``"cupy"`` is selected, ensure that CuPy is installed and a compatible
      GPU is available.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(4)
    >>> L, w = get_graph_lspectrum(G, library='numpy')
    >>> print(L)
    [[ 1 -1  0  0]
     [-1  2 -1  0]
     [ 0 -1  2 -1]
     [ 0  0 -1  1]]
    >>> print(w)
    [0. 0.58578644 2. 3.41421356]
    """
    L = nx.laplacian_matrix(G)

    match library:
        case 'numpy':
            L = L.toarray()
            w = np.linalg.eigvals(L)
        case 'networkx'|'nx'|'scipy'|'sp'|'cupy'|'cp':
            match library:
                case 'networkx'|'nx':
                    w = nx.laplacian_spectrum(G)
                case 'scipy'|'sp':
                    from scipy.linalg import eigvals
                    w = eigvals(L)
                case 'cupy':
                    import cupy as cp
                    L_gpu = cp.asarray(L)
                    w = cp.linalg.eigvals(L_gpu).get()
            L = L.toarray()
        case _:
            raise ValueError("Unsupported library. Choose from 'numpy', 'networkx', 'scipy', or 'cupy'.")

    return L, w


def get_graph_lspectrum_rw(
    G: nx.Graph,
    is_signed: bool = False,
) -> tuple[NDArray, NDArray]:
    """
    Compute the random-walk normalized Laplacian and its spectrum.

    Parameters
    ----------
    G : nx.Graph
        The input graph.
    is_signed : bool, default False
        If True, use the signed adjacency to build the normalized Laplacian
        and compute eigenvalues directly from it.

    Returns
    -------
    tuple[NDArray, NDArray]
        A tuple containing the normalized Laplacian and its eigenvalues.

    Notes
    -----
    When ``is_signed=False``, eigenvalues are computed with
    ``networkx.laplacian_spectrum`` (combinatorial Laplacian), which may not
    match the normalized Laplacian returned as ``L``.
    """
    A = nx.adjacency_matrix(G).toarray()
    D = np.diag(np.abs(A).sum(axis=1))
    L = np.eye(D.shape[0]) - fractional_matrix_power(
        D, -0.5
    ) @ A @ fractional_matrix_power(D, -0.5)
    if is_signed:
        w = np.linalg.eigvals(L)
    else:
        w = nx.laplacian_spectrum(G)
    return L, w


def compute_laplacian_properties(
    G: nx.Graph,
    tau: Optional[float] = None,
) -> tuple[NDArray, NDArray, NDArray, NDArray, float]:
    """
    Computes the Laplacian spectrum, Laplacian matrix, rho, and Trho for graph G.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    tau : float, optional
        Scale parameter for the matrix exponential. If None, uses
        ``1 / max(spectrum)``.

    Returns
    -------
    spectrum : ndarray
        The Laplacian spectrum of G.
    L : ndarray
        The dense Laplacian matrix of G.
    rho : ndarray
        Matrix computed as ``expm(-tau * L)`` normalized by its trace.
    Trho : ndarray
        Symmetric matrix obtained from the element-wise inverse of ``rho``,
        with its diagonal set to 0.
    tau : float
        The value of ``tau`` used in the computation.
    """
    spectrum = nx.laplacian_spectrum(G)
    L_sparse = nx.laplacian_matrix(G)
    L = L_sparse.todense()
    tau = tau or 1/max(spectrum)
    num = expm(-tau * L)
    den = np.trace(num)
    rho = num / den
    Trho = np.copy(1.0 / rho)
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0)
    
    return spectrum, L, rho, Trho, tau

