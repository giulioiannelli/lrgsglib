"""Spectral utilities for signed graphs."""

from __future__ import annotations

import numpy as np
from networkx import Graph, NetworkXError, to_numpy_array, to_scipy_sparse_array
from networkx.drawing.layout import _process_params, rescale_layout
from networkx.utils.backends import _dispatchable
from scipy.sparse import csr_array, spdiags
from scipy.sparse.linalg import eigsh
from typing import List

__all__ = [
    "signed_laplacian_matrix",
    "signed_spectral_layout",
]


@_dispatchable(edge_attrs="weight")
def signed_laplacian_matrix(
    G: Graph, nodelist: List | None = None, weight: str = "weight"
) -> csr_array:
    """Return the signed Laplacian matrix of ``G``."""
    nodelist = nodelist or list(G)
    adj = to_scipy_sparse_array(G, nodelist=nodelist, weight=weight, format="csr")
    deg = csr_array(spdiags(abs(adj).sum(axis=1), 0, *adj.shape, format="csr"))
    return deg - adj


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
