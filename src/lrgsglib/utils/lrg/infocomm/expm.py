"""Entropy computation using expm_multiply.

This module provides entropy computation using scipy's expm_multiply
for efficient matrix exponential computation on sparse matrices.
"""

from typing import Optional, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_array, csr_matrix
from scipy.sparse.linalg import expm_multiply

from ...tools.loglib import get_logger
from ._common import _normalize_entropy_profile

__all__ = [
    "compute_entropy_observables_expm_multiply",
]

logger = get_logger(__name__)


def compute_entropy_observables_expm_multiply(
    L: Union[np.ndarray, csr_matrix, csr_array],
    num_nodes: int,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    num_samples: int = 30,
    seed: Optional[int] = None,
    typf: type = np.float64,
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
    verbose: bool = False,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Compute entropy observables using expm_multiply (bypasses eigenvalues).

    Uses Hutchinson stochastic trace estimator with Rademacher probes to
    compute von Neumann entropy without explicit eigenvalue decomposition.
    Suitable for large sparse Laplacians where eigendecomposition is expensive.

    Mathematical Formula
    --------------------
    S(τ) = log Z(τ) + τ · Tr[L exp(-τL)] / Z(τ)

    where:
    - Z(τ) = Tr[exp(-τL)] (partition function)
    - Both traces estimated via Hutchinson estimator with Rademacher probes

    Parameters
    ----------
    L : sparse matrix or ndarray
        Signed Laplacian matrix (must be positive semi-definite).
    num_nodes : int
        Number of nodes (for normalization).
    steps : int, optional
        Number of logarithmic time samples. Default: 600.
    t1, t2 : int, optional
        Time range: τ ∈ [10^t1, 10^t2]. Default: [-2, 5].
    num_samples : int, optional
        Number of Rademacher probes for Hutchinson estimator. Default: 30.
        Error decreases as O(1/√num_samples).
    seed : int, optional
        Random seed for reproducibility of probe vectors.
    typf : type, optional
        Float precision. Default: np.float64.
    entropy_norm : str, optional
        Normalization mode: "standard" (raw) or "complement" (1 - S). Default: "complement".
    specific_heat_scale : str, optional
        Scaling: "logN" (remove size dependence) or "none". Default: "logN".
    verbose : bool, optional
        Deprecated. Use ``lrgsglib.utils.tools.loglib.enable_logging()`` instead.
        Default: False.

    Returns
    -------
    normalized_entropy : NDArray
        Entropy profile over time grid. Shape: (steps,).
    specific_heat : NDArray
        Entropy derivative (generalized specific heat). Shape: (steps,).
    variance_profile : NDArray
        Variance profile (zeros in this mode). Shape: (steps,).
    time_grid : NDArray
        Logarithmic time grid. Shape: (steps,).

    Notes
    -----
    - Requires scipy >= 1.6.0 for expm_multiply.
    - Variance computation not implemented (would require Tr[L² exp(-τL)]).
    - For positive semi-definite Laplacians, exp(-τL) decays exponentially.
    - Same probe vectors used across all τ to reduce noise in entropy derivative.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
    >>> G = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, p3=0.6, p4=0.8,
    ...                                 fraction=0.4, iterations=5)
    >>> G.upd_graph_matrices()  # Ensure slp is computed
    >>> S, C, V, tau = compute_entropy_observables_expm_multiply(
    ...     L=G.slp, num_nodes=G.N, steps=200, num_samples=50, seed=42
    ... )
    >>> print(f"Entropy at τ=1: {S[np.argmin(np.abs(tau - 1.0))]:.4f}")

    References
    ----------
    Hutchinson, M. F. (1990). A stochastic estimator of the trace of the
    influence matrix for Laplacian smoothing splines. Communications in
    Statistics - Simulation and Computation, 19(2), 433-450.
    """
    import scipy.sparse

    if steps < 1:
        raise ValueError("steps must be at least 1")

    # Generate time grid
    time_grid = np.logspace(t1, t2, steps, dtype=typf)

    # Ensure L is sparse (expm_multiply works best with sparse matrices)
    if not scipy.sparse.issparse(L):
        L = scipy.sparse.csr_matrix(L, dtype=typf)
    else:
        L = L.astype(typf, copy=False)

    N = L.shape[0]
    if N != num_nodes:
        raise ValueError(f"L shape {L.shape} incompatible with num_nodes={num_nodes}")

    # Generate Rademacher probes (±1 random vectors)
    rng = np.random.default_rng(seed)
    Z = rng.integers(0, 2, size=(N, num_samples), dtype=np.int8)
    Z = (2 * Z - 1).astype(typf, copy=False)  # Convert {0,1} -> {-1,+1}

    # Allocate output array
    entropy_profile = np.zeros(steps, dtype=typf)

    # Main loop: compute entropy for each time point
    tau_cutoff = 5000  # Aggressive cutoff for tractable computation
    _logged_cutoff = False

    for idx, tau in enumerate(time_grid):
        # For very large tau, skip expensive computation
        if tau > tau_cutoff:
            entropy_profile[idx] = 0.0
            if not _logged_cutoff:
                logger.debug(
                    "Using limiting value S=0 for tau > %d (avoids slow computation)",
                    tau_cutoff,
                )
                _logged_cutoff = True
            continue

        # Matrix exponential action: Y = exp(-τL) Z
        A = (-tau) * L
        try:
            Y = expm_multiply(A, Z)  # Shape (N, num_samples)
        except Exception as e:
            logger.warning("expm_multiply failed at tau=%.2e: %s", tau, e)
            entropy_profile[idx] = np.nan
            continue

        # Hutchinson trace estimates
        # Z(τ) = Tr[exp(-τL)] ≈ (1/m) Σ z_j^T y_j
        trK = np.mean(np.sum(Z * Y, axis=0))

        # Tr[L exp(-τL)] ≈ (1/m) Σ z_j^T (L y_j)
        LY = L @ Y
        trLK = np.mean(np.sum(Z * LY, axis=0))

        # Sanity check
        if not np.isfinite(trK) or trK <= 0:
            logger.warning("Invalid partition function Z(tau=%.2e) = %.2e", tau, trK)
            entropy_profile[idx] = np.nan
            continue

        # Compute entropy: S(τ) = log Z(τ) + τ · Tr[L exp(-τL)] / Z(τ)
        entropy_profile[idx] = np.log(trK) + tau * (trLK / trK)

    # Normalize entropy
    log_N = np.log(typf(num_nodes)) if num_nodes > 0 else typf(1)
    entropy_profile = entropy_profile / log_N
    normalized_entropy = _normalize_entropy_profile(entropy_profile, entropy_norm)

    # Compute specific heat (entropy derivative with respect to log τ)
    specific_heat = np.gradient(normalized_entropy, np.log(time_grid))
    scale = specific_heat_scale.lower()
    if scale == "logn":
        specific_heat = log_N * specific_heat
    elif scale != "none":
        raise ValueError("specific_heat_scale must be 'logN' or 'none'")

    # Variance not computed in this mode
    variance_profile = np.zeros_like(time_grid, dtype=typf)

    return normalized_entropy, specific_heat, variance_profile, time_grid
