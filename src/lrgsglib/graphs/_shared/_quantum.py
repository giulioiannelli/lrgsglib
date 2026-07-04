"""Engine-agnostic quantum propagator mixin functions."""

from typing import Union
import numpy as np
from numpy.typing import NDArray


def _eigV_column_major(self) -> NDArray:
    """Return eigV in column-major format (each column = one eigenvector).

    NX stores eigV row-major (transposed) by default; GT stores column-major.
    The quantum utility functions expect column-major.
    """
    eigV = self.eigV
    if getattr(self, "_eigV_is_transposed", False):
        return eigV.T
    return eigV


def compute_quantum_propagator(
    self,
    t: float,
    full_matrix: bool = False,
) -> Union[tuple[NDArray, NDArray], NDArray]:
    """Compute quantum propagator U(t) = exp(-i*t*L) for this graph.

    Requires eigendecomposition (``self.eigv``, ``self.eigV``).
    """
    if not hasattr(self, "eigv") or not hasattr(self, "eigV"):
        raise ValueError(
            "Eigendecomposition required. "
            "Call compute_laplacian_spectrum_weigV() first."
        )
    from ...utils.lrg.quantum import compute_quantum_propagator_spectral

    return compute_quantum_propagator_spectral(
        self.eigv, _eigV_column_major(self), t, full_matrix=full_matrix
    )


def quantum_walk_probabilities(
    self,
    t: float,
    init_node: int,
) -> NDArray:
    """Compute quantum walk probability distribution at time t.

    Returns ``P_j(t) = |<j|U(t)|i>|^2`` starting from ``init_node``.
    """
    if not hasattr(self, "eigv") or not hasattr(self, "eigV"):
        raise ValueError(
            "Eigendecomposition required. "
            "Call compute_laplacian_spectrum_weigV() first."
        )
    if init_node < 0 or init_node >= self.N:
        raise ValueError(f"init_node={init_node} out of range [0, {self.N})")

    from ...utils.lrg.quantum import quantum_probability_distribution

    return quantum_probability_distribution(
        self.eigv, _eigV_column_major(self), t, init_node
    )


def quantum_observables_time_series(
    self,
    t_array: NDArray,
    init_type: str = "localized",
    init_node: int = 0,
) -> dict[str, NDArray]:
    """Compute quantum observables (entropy, coherence, purity) over time."""
    if not hasattr(self, "eigv") or not hasattr(self, "eigV"):
        raise ValueError(
            "Eigendecomposition required. "
            "Call compute_laplacian_spectrum_weigV() first."
        )
    from ...utils.lrg.quantum import compute_quantum_observables_from_eigenvalues

    return compute_quantum_observables_from_eigenvalues(
        self.eigv, _eigV_column_major(self), self.N,
        t_array, init_type, init_node
    )
