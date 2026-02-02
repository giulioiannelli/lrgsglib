"""Information-theoretic entropy observables for Laplacian spectra.

This package provides multiple approaches for computing entropy-related
observables from graph Laplacian spectra:

- **classical**: Exact computation from eigenvalues (recommended for small graphs)
- **slq**: Stochastic Lanczos Quadrature (for large matrices without eigendecomposition)
- **expm**: Matrix exponential approach using expm_multiply
- **renyi**: Generalized Rényi and Tsallis entropies
- **distances**: Spectral distances and ultrametric matrices
"""

# Core classical approach
from .classical import (
    entropy,
    compute_entropy_observables_from_eigenvalues,
)

# Alternative computation methods
from .slq import compute_entropy_observables_slq
from .expm import compute_entropy_observables_expm_multiply

# Generalized entropies
from .renyi import compute_renyi_observables_from_eigenvalues

# Distance functions
from .distances import (
    extract_ultrametric_matrix,
    lapl_dists,
)

__all__ = [
    # Classical
    "entropy",
    "compute_entropy_observables_from_eigenvalues",
    # SLQ
    "compute_entropy_observables_slq",
    # expm
    "compute_entropy_observables_expm_multiply",
    # Renyi
    "compute_renyi_observables_from_eigenvalues",
    # Distances
    "extract_ultrametric_matrix",
    "lapl_dists",
]
