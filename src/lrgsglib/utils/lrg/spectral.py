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
from ...config.const import (
    SG_LAPL_SIGNED,
    SG_LAPL_RW,
    SG_LAPL_SYM,
    SG_LAPL_TYPES,
    SG_LAPL_DEFAULT_TYPE,
    SG_LAPL_RW_IMAG_TOL,
)
#
__all__ = [
    "get_graph_lspectrum",
    "get_graph_lspectrum_rw",
    "compute_laplacian_properties",
]


def get_graph_lspectrum(
    G: nx.Graph,
    library: str = "numpy",
    signed: bool = False,
    laplacian_type: str = SG_LAPL_DEFAULT_TYPE,
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
    signed : bool, optional
        If True, compute the signed Laplacian (L = D - A where D uses absolute
        degrees and A preserves edge signs). This is appropriate for graphs with
        negative edge weights representing antiferromagnetic/repulsive interactions.
        If False (default), uses NetworkX's standard unsigned Laplacian.
        Only consulted when ``laplacian_type == 'signed'``.
    laplacian_type : str, optional
        Which signed Laplacian to build (Kunegis 2010):
        ``'signed'`` (default) -> combinatorial (honours the ``signed`` flag);
        ``'sym'`` -> symmetric normalized ``L_sym = I - D_s^-1/2 A D_s^-1/2``;
        ``'rw'`` -> random-walk ``L_rw = I - D_s^-1 A`` (non-symmetric; its
        eigenvalues are real and equal those of ``'sym'``).

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
    - When ``signed=False``, the Laplacian is computed using NetworkX's
      ``laplacian_matrix``, which ignores edge signs.
    - When ``signed=True``, uses ``signed_laplacian_matrix`` from
      ``graphs.nx.funcs.spectral``, which correctly handles negative edge weights.
    - For signed graphs (with negative edges), the signed Laplacian may have
      negative eigenvalues, unlike the standard Laplacian which is always PSD.
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

    For signed graphs with negative edges:

    >>> G = nx.Graph()
    >>> G.add_edge(0, 1, weight=-1)  # antiferromagnetic edge
    >>> G.add_edge(1, 2, weight=1)
    >>> L, w = get_graph_lspectrum(G, library='numpy', signed=True)
    """
    if laplacian_type in (SG_LAPL_SYM, SG_LAPL_RW):
        from lrgsglib.graphs.nx.funcs.spectral import (
            signed_sym_laplacian_matrix,
            signed_rw_laplacian_matrix,
        )

        is_sym = laplacian_type == SG_LAPL_SYM
        builder = signed_sym_laplacian_matrix if is_sym else signed_rw_laplacian_matrix
        L = builder(G).toarray()
        match library:
            case "scipy" | "sp":
                from scipy.linalg import eigvalsh as _eigvalsh, eigvals as _eigvals
                w = _eigvalsh(L) if is_sym else _eigvals(L)
            case "cupy" | "cp":
                import cupy as cp

                Lg = cp.asarray(L)
                w = (cp.linalg.eigvalsh(Lg) if is_sym else cp.linalg.eigvals(Lg)).get()
            case _:
                # numpy / networkx / default: nx has no normalized signed
                # spectrum, so compute directly from L.
                w = np.linalg.eigvalsh(L) if is_sym else np.linalg.eigvals(L)
        if not is_sym:
            # L_rw is isospectral to the symmetric L_sym -> eigenvalues are real
            imag = float(np.max(np.abs(np.imag(w)))) if np.iscomplexobj(w) else 0.0
            if imag > SG_LAPL_RW_IMAG_TOL:
                import warnings

                warnings.warn(
                    f"rw Laplacian eigenvalues have |Im|={imag:.2e} > "
                    f"{SG_LAPL_RW_IMAG_TOL:.1e}; casting to real.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            w = np.real(w)
        return L, w

    if laplacian_type != SG_LAPL_SIGNED:
        raise ValueError(
            f"laplacian_type must be one of {SG_LAPL_TYPES}, got {laplacian_type!r}"
        )

    if signed:
        # Import signed Laplacian from canonical location
        from lrgsglib.graphs.nx.funcs.spectral import signed_laplacian_matrix

        L = signed_laplacian_matrix(G)
    else:
        L = nx.laplacian_matrix(G)

    match library:
        case "numpy":
            L = L.toarray()
            w = np.linalg.eigvals(L)
        case "networkx" | "nx" | "scipy" | "sp" | "cupy" | "cp":
            match library:
                case "networkx" | "nx":
                    if signed:
                        # networkx.laplacian_spectrum ignores signs, compute directly
                        L = L.toarray()
                        w = np.linalg.eigvals(L)
                    else:
                        w = nx.laplacian_spectrum(G)
                        L = L.toarray()
                case "scipy" | "sp":
                    from scipy.linalg import eigvals

                    L = L.toarray()
                    w = eigvals(L)
                case "cupy" | "cp":
                    import cupy as cp

                    L = L.toarray()
                    L_gpu = cp.asarray(L)
                    w = cp.linalg.eigvals(L_gpu).get()
        case _:
            raise ValueError(
                "Unsupported library. Choose from 'numpy', 'networkx', 'scipy', or 'cupy'."
            )

    return L, w


def get_graph_lspectrum_rw(
    G: nx.Graph,
    is_signed: bool = False,
) -> tuple[NDArray, NDArray]:
    """
    Symmetric normalized signed Laplacian and its spectrum (thin alias).

    .. deprecated::
        Prefer :func:`get_graph_lspectrum` with ``laplacian_type='sym'`` (this
        matrix) or ``laplacian_type='rw'`` (the true non-symmetric random-walk
        operator). Retained for backward compatibility.

    Returns ``L_sym = I - D_s^{-1/2} A D_s^{-1/2}`` (Kunegis 2010, §3.4) with the
    signed degree ``D_s`` and its real eigenvalues. The ``is_signed`` flag is
    accepted for signature compatibility but no longer changes the result: the
    spectrum is always computed from the returned (signed-degree) matrix.
    """
    from lrgsglib.graphs.nx.funcs.spectral import signed_sym_laplacian_matrix

    L = signed_sym_laplacian_matrix(G).toarray()
    w = np.linalg.eigvalsh(L)
    return L, w


def compute_laplacian_properties(
    G: nx.Graph,
    tau: Optional[float] = None,
    signed: bool = False,
) -> tuple[NDArray, NDArray, NDArray, NDArray, float]:
    """
    Computes the Laplacian spectrum, Laplacian matrix, rho, and Trho for graph G.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    tau : float, optional
        Scale parameter for the matrix exponential. If None, uses
        ``1 / max(abs(spectrum))``.
    signed : bool, optional
        If True, compute the signed Laplacian (L = D - A where D uses absolute
        degrees and A preserves edge signs). Default is False.

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

    Notes
    -----
    When ``signed=True``, the signed Laplacian may have negative eigenvalues.
    The tau parameter is computed from ``max(abs(spectrum))`` to handle this case.
    """
    if signed:
        from lrgsglib.graphs.nx.funcs.spectral import signed_laplacian_matrix

        L_sparse = signed_laplacian_matrix(G)
        L = np.asarray(L_sparse.todense())
        spectrum = np.linalg.eigvalsh(L)
    else:
        spectrum = nx.laplacian_spectrum(G)
        L_sparse = nx.laplacian_matrix(G)
        L = np.asarray(L_sparse.todense())

    # Use absolute value of spectrum for tau calculation (handles negative eigenvalues)
    tau = tau or 1 / max(np.abs(spectrum))
    num = expm(-tau * L)
    den = np.trace(num)
    rho = num / den
    Trho = np.copy(1.0 / rho)
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0)

    return spectrum, L, rho, Trho, tau

