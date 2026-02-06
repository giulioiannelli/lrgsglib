"""Spectral computation methods for GraphOfGraphs structures.

These methods exploit the product structure of graph-of-graphs to compute
eigenvalues efficiently by separating base and fiber spectra.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from .GraphOfGraphsNX import GraphOfGraphsNX


def compute_base_spectrum(self: "GraphOfGraphsNX", typf: type = np.float64) -> np.ndarray:
    """Compute eigenvalues of the base graph Laplacian.

    Parameters
    ----------
    typf : type
        Data type for computation

    Returns
    -------
    np.ndarray
        Base graph eigenvalues in ascending order
    """
    if self._base_eigenvalues is not None:
        return self._base_eigenvalues

    base_lap = nx.laplacian_matrix(self._base_graph_instance.G).astype(typf).todense()
    self._base_eigenvalues = np.linalg.eigvalsh(base_lap)
    return self._base_eigenvalues


def compute_fiber_laplacian(self: "GraphOfGraphsNX", typf: type = np.float64) -> np.ndarray:
    """Compute the fiber graph Laplacian matrix.

    Parameters
    ----------
    typf : type
        Data type for computation

    Returns
    -------
    np.ndarray
        Dense Laplacian matrix of the fiber graph
    """
    if self._fiber_laplacian is not None:
        return self._fiber_laplacian

    fiber_lap = nx.laplacian_matrix(self._fiber_graph_instance.G).astype(typf).todense()
    self._fiber_laplacian = np.asarray(fiber_lap, dtype=typf)
    return self._fiber_laplacian


def compute_separated_spectrum(
    self: "GraphOfGraphsNX",
    typf: type = np.float64,
    backend: str = "numpy",
) -> np.ndarray:
    """Compute the full spectrum using separated computation.

    This method uses the product structure of graph-of-graphs to compute the
    spectrum more efficiently than diagonalizing the full Laplacian.
    Complexity: O(N_base³ + N_base * N_fiber³) vs O((N_base * N_fiber)³).

    Only valid when can_use_separated_spectrum() returns True.

    Parameters
    ----------
    typf : type
        Data type for computation
    backend : str
        Backend for eigenvalue computation (currently only 'numpy')

    Returns
    -------
    np.ndarray
        Full sorted spectrum of the structure

    Raises
    ------
    ValueError
        If separated spectrum computation is not valid for this structure

    Notes
    -----
    The eigenvalues are stored in self.eigv for compatibility with
    SignedGraph spectral methods.
    """
    if not self.can_use_separated_spectrum():
        raise ValueError(
            "Separated spectrum not valid for heterogeneous fibers or "
            "non-uniform anchor indices. Use standard spectral methods."
        )

    base_eigenvalues = compute_base_spectrum(self, typf=typf)
    fiber_lap = compute_fiber_laplacian(self, typf=typf)

    # Get anchor index (uniform for all fibers when separated spectrum is valid)
    anchor_idx = self._anchor_indices[0]

    all_eigenvalues = []
    for lambda_base in base_eigenvalues:
        # Create modified fiber Laplacian: L_fiber + λ_base * e_anchor * e_anchor^T
        impure_lap = fiber_lap.copy()
        impure_lap[anchor_idx, anchor_idx] += lambda_base
        impure_eigenvalues = np.linalg.eigvalsh(impure_lap)
        all_eigenvalues.extend(impure_eigenvalues)

    eigenvalues = np.sort(np.array(all_eigenvalues, dtype=typf))
    self.eigv = eigenvalues  # Store for SignedGraph compatibility
    return eigenvalues


def get_base_fiber_dimensions(self: "GraphOfGraphsNX") -> tuple[int, int]:
    """Get the dimensions of base and fiber components.

    Returns
    -------
    tuple[int, int]
        (N_base, N_fiber)
    """
    return (self._n_base, self._n_fiber)
