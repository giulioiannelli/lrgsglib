"""Spectral utilities for signed graphs."""

from __future__ import annotations

import numpy as np
from networkx import Graph, NetworkXError, to_numpy_array, to_scipy_sparse_array
from networkx.drawing.layout import _process_params, rescale_layout
from networkx.utils.backends import _dispatchable, _registered_algorithms
from scipy.sparse import csr_array, spdiags, identity as sp_identity
from scipy.sparse.linalg import eigsh
from typing import List

__all__ = [
    "signed_laplacian_matrix",
    "signed_sym_laplacian_matrix",
    "signed_rw_laplacian_matrix",
    "signed_spectral_layout",
]


def _signed_adj_and_absdeg(G, nodelist=None, weight="weight"):
    """Signed adjacency ``A`` (CSR) and absolute-degree vector ``D_s[i,i]``.

    ``D_s[i] = sum_j |A_ij|`` is the signed degree of Kunegis 2010 (Eq. 3.8),
    shared by the combinatorial, symmetric-normalized and random-walk signed
    Laplacians.
    """
    nodelist = nodelist or list(G)
    adj = to_scipy_sparse_array(G, nodelist=nodelist, weight=weight, format="csr")
    deg = np.asarray(abs(adj).sum(axis=1)).ravel()
    return adj, deg

if "signed_laplacian_matrix" in _registered_algorithms:
    def _dispatchable_signed(func):
        return func
else:
    _dispatchable_signed = _dispatchable(edge_attrs="weight")


@_dispatchable_signed
def signed_laplacian_matrix(
    G: Graph, nodelist: List | None = None, weight: str = "weight"
) -> csr_array:
    """Return the signed Laplacian matrix of ``G``."""
    nodelist = nodelist or list(G)
    adj = to_scipy_sparse_array(G, nodelist=nodelist, weight=weight, format="csr")
    deg = csr_array(spdiags(abs(adj).sum(axis=1), 0, *adj.shape, format="csr"))
    return deg - adj


def signed_sym_laplacian_matrix(
    G: Graph, nodelist: List | None = None, weight: str = "weight"
) -> csr_array:
    """Return the symmetric normalized signed Laplacian of ``G``.

    ``L_sym = I - D_s^{-1/2} A D_s^{-1/2}`` (Kunegis 2010, §3.4). Symmetric and
    positive-semidefinite, spectrum in ``[0, 2]``. Isolated nodes (zero signed
    degree) get a zero inverse, leaving that row equal to the identity row
    (eigenvalue 1).
    """
    adj, deg = _signed_adj_and_absdeg(G, nodelist=nodelist, weight=weight)
    n = adj.shape[0]
    inv_sqrt = np.zeros_like(deg, dtype=float)
    nz = deg > 0
    inv_sqrt[nz] = 1.0 / np.sqrt(deg[nz])
    dinv_sqrt = csr_array(spdiags(inv_sqrt, 0, n, n, format="csr"))
    return csr_array(sp_identity(n, format="csr")) - dinv_sqrt @ adj @ dinv_sqrt


def signed_rw_laplacian_matrix(
    G: Graph, nodelist: List | None = None, weight: str = "weight"
) -> csr_array:
    """Return the random-walk normalized signed Laplacian of ``G``.

    ``L_rw = I - D_s^{-1} A`` (Kunegis 2010, §3.3). This matrix is NOT symmetric
    (use a general eigensolver), but it is similar to ``L_sym`` via
    ``L_rw = D_s^{-1/2} L_sym D_s^{1/2}`` and hence isospectral: its eigenvalues
    are real and lie in ``[0, 2]``. Isolated nodes (zero signed degree) get a
    zero inverse, leaving that row equal to the identity row (eigenvalue 1).
    """
    adj, deg = _signed_adj_and_absdeg(G, nodelist=nodelist, weight=weight)
    n = adj.shape[0]
    inv = np.zeros_like(deg, dtype=float)
    nz = deg > 0
    inv[nz] = 1.0 / deg[nz]
    dinv = csr_array(spdiags(inv, 0, n, n, format="csr"))
    return csr_array(sp_identity(n, format="csr")) - dinv @ adj


def signed_spectral_layout(
    G: Graph, weight: str = "weight", scale: float = 1, center=None, dim: int = 2
):
    """Position nodes using eigenvectors of the signed Laplacian."""
    G, center = _process_params(G, center, dim)

    if len(G) <= 2:
        if len(G) == 0:
            pos = np.array([])
        elif len(G) == 1:
            pos = np.array([center])
        else:
            pos = np.array([np.zeros(dim), np.array(center) * 2.0])
        return dict(zip(G, pos))
    try:
        if len(G) < 500:
            raise ValueError
        A = to_scipy_sparse_array(G, weight=weight, dtype="d")
        if G.is_directed():
            A = A + A.transpose()
        pos = _sparse_spectral_signed(A, dim)
    except (ImportError, ValueError):
        A = to_numpy_array(G, weight=weight)
        if G.is_directed():
            A += A.transpose()
        pos = _spectral_signed(A, dim)

    pos = rescale_layout(pos, scale=scale) + center
    return dict(zip(G, pos))


def _spectral_signed(A, dim=2):
    try:
        nnodes, _ = A.shape
    except AttributeError as err:
        msg = "spectral() takes an adjacency matrix as input"
        raise NetworkXError(msg) from err

    D = np.identity(nnodes, dtype=A.dtype) * abs(A).sum(axis=1)
    L = D - A

    eigenvalues, eigenvectors = np.linalg.eig(L)
    index = np.argsort(eigenvalues)[1 : dim + 1]
    return np.real(eigenvectors[:, index])


def _sparse_spectral_signed(A, dim=2):
    try:
        nnodes, _ = A.shape
    except AttributeError as err:
        msg = "sparse_spectral() takes an adjacency matrix as input"
        raise NetworkXError(msg) from err

    D = csr_array(spdiags(abs(A).sum(axis=1), 0, nnodes, nnodes))
    L = D - A

    k = dim + 1
    ncv = max(2 * k + 1, int(np.sqrt(nnodes)))
    eigenvalues, eigenvectors = eigsh(L, k, which="SM", ncv=ncv)
    index = np.argsort(eigenvalues)[1:k]
    return np.real(eigenvectors[:, index])
