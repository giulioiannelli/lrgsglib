"""Spectral computation methods for Dirac lattice structures.

These methods exploit the product structure of Dirac graphs to compute
eigenvalues efficiently by separating base and fiber spectra.
"""

import networkx as nx
import numpy as np


def compute_base_spectrum(self, typf: type = np.float64) -> np.ndarray:
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
    if self.base_eigenvalues is not None:
        return self.base_eigenvalues

    base_lap = nx.laplacian_matrix(self.base_graph).astype(typf).todense()
    self.base_eigenvalues = np.linalg.eigvalsh(base_lap)
    return self.base_eigenvalues


def compute_fiber_laplacian(self, typf: type = np.float64) -> np.ndarray:
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
    if self.fiber_laplacian is not None:
        return self.fiber_laplacian

    fiber_lap = nx.laplacian_matrix(self.fiber_graph).astype(typf).todense()
    self.fiber_laplacian = np.asarray(fiber_lap, dtype=typf)
    return self.fiber_laplacian


def compute_dirac_spectrum_separated(
    self, typf: type = np.float64, backend: str = 'numpy'
) -> np.ndarray:
    """Compute the full spectrum of a Dirac structure efficiently.

    This method uses the product structure of Dirac graphs to compute the
    spectrum more efficiently than diagonalizing the full Laplacian.
    Complexity: O(N_base^3 + N_base * N_fiber^3) vs O((N_base * N_fiber)^3).

    Parameters
    ----------
    typf : type
        Data type for computation
    backend : str
        Backend for eigenvalue computation (currently only 'numpy')

    Returns
    -------
    np.ndarray
        Full sorted spectrum of the Dirac structure

    Notes
    -----
    The eigenvalues are stored in self.eigv for compatibility with
    SignedGraph spectral methods.
    """
    base_eigenvalues = self.compute_base_spectrum(typf=typf)
    fiber_lap = self.compute_fiber_laplacian(typf=typf)

    all_eigenvalues = []
    for lambda_base in base_eigenvalues:
        # Create modified fiber Laplacian: L_fiber + lambda_base * e1e1^T
        impure_lap = fiber_lap.copy()
        impure_lap[0, 0] += lambda_base
        impure_eigenvalues = np.linalg.eigvalsh(impure_lap)
        all_eigenvalues.extend(impure_eigenvalues)

    eigenvalues = np.sort(np.array(all_eigenvalues, dtype=typf))
    self.eigv = eigenvalues  # Store for SignedGraph compatibility
    return eigenvalues


def get_base_fiber_dimensions(self) -> tuple[int, int]:
    """Get the dimensions of base and fiber components.

    Returns
    -------
    tuple[int, int]
        (base_nodes, fiber_nodes)
    """
    return (
        self.dirac_structure['base_nodes'],
        self.dirac_structure['fiber_nodes']
    )
