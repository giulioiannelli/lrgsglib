"""Stochastic Lanczos Quadrature (SLQ) entropy computation.

This module provides entropy computation using SLQ for large matrices
where full eigenvalue decomposition is impractical.
"""

from typing import Callable, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh_tridiagonal

from ._common import _normalize_entropy_profile

__all__ = [
    "compute_entropy_observables_slq",
]


def _lanczos_tridiagonal(
    matvec: Callable[[NDArray], NDArray],
    n: int,
    m: int,
    v0: NDArray,
    typf: type,
) -> Tuple[NDArray, NDArray]:
    """Lanczos algorithm to build tridiagonal matrix.

    Parameters
    ----------
    matvec : Callable
        Matrix-vector product function.
    n : int
        Matrix dimension.
    m : int
        Number of Lanczos steps.
    v0 : NDArray
        Starting vector.
    typf : type
        Floating-point type.

    Returns
    -------
    Tuple[NDArray, NDArray]
        Diagonal (alphas) and off-diagonal (betas) of tridiagonal matrix.
    """
    v_prev = np.zeros(n, dtype=typf)
    v = v0 / np.linalg.norm(v0)
    alphas = np.zeros(m, dtype=typf)
    betas = np.zeros(m - 1, dtype=typf)
    beta = typf(0)

    for j in range(m):
        w = matvec(v)
        if j > 0:
            w = w - beta * v_prev
        alpha = np.dot(v, w)
        w = w - alpha * v
        beta = np.linalg.norm(w)

        alphas[j] = alpha
        if j < m - 1:
            betas[j] = beta

        if beta == 0:
            return alphas[: j + 1], betas[:j]

        v_prev, v = v, w / beta

    return alphas, betas


def _slq_trace_estimates(
    laplacian: NDArray,
    time_grid: NDArray,
    num_samples: int,
    lanczos_steps: int,
    seed: Optional[int],
    typf: type,
) -> Tuple[NDArray, NDArray, NDArray]:
    """Estimate traces via stochastic Lanczos quadrature.

    Parameters
    ----------
    laplacian : NDArray
        Laplacian matrix.
    time_grid : NDArray
        Time points for evaluation.
    num_samples : int
        Number of random probe vectors.
    lanczos_steps : int
        Number of Lanczos iterations.
    seed : int, optional
        Random seed.
    typf : type
        Floating-point type.

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray]
        Estimates of Tr[exp(-τL)], Tr[L exp(-τL)], Tr[L² exp(-τL)].
    """
    n = laplacian.shape[0]
    rng = np.random.default_rng(seed)
    trace_exp = np.zeros_like(time_grid, dtype=typf)
    trace_l_exp = np.zeros_like(time_grid, dtype=typf)
    trace_l2_exp = np.zeros_like(time_grid, dtype=typf)

    matvec = laplacian.dot if hasattr(laplacian, "dot") else lambda x: laplacian @ x

    for _ in range(num_samples):
        r = rng.choice([-1.0, 1.0], size=n).astype(typf)
        norm_sq = np.dot(r, r)
        alphas, betas = _lanczos_tridiagonal(
            matvec=matvec,
            n=n,
            m=lanczos_steps,
            v0=r,
            typf=typf,
        )
        evals, evecs = eigh_tridiagonal(alphas, betas)
        weights = evecs[0, :] ** 2

        exp_terms = np.exp(-np.outer(time_grid, evals))
        w0 = exp_terms @ weights
        w1 = exp_terms @ (weights * evals)
        w2 = exp_terms @ (weights * (evals**2))

        trace_exp += norm_sq * w0
        trace_l_exp += norm_sq * w1
        trace_l2_exp += norm_sq * w2

    inv_samples = typf(1.0 / num_samples)
    return trace_exp * inv_samples, trace_l_exp * inv_samples, trace_l2_exp * inv_samples


def compute_entropy_observables_slq(
    laplacian: NDArray,
    num_nodes: int,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    typf: type = np.float64,
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
    lanczos_steps: int = 50,
    num_samples: int = 30,
    seed: Optional[int] = None,
    tau_cutoff: Optional[float] = None,
    auto_cutoff_factor: float = 10.0,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Approximate entropy observables with stochastic Lanczos quadrature.

    This avoids explicit eigenvalue computation by estimating traces of
    exp(-tau L), L exp(-tau L), and L^2 exp(-tau L) with Hutchinson sampling.

    For large tau where SLQ becomes numerically unstable, uses analytical
    limiting values. The cutoff is either user-specified or automatically
    estimated from the spectral gap.

    Parameters
    ----------
    laplacian : NDArray
        Laplacian matrix (dense or sparse).
    num_nodes : int
        Number of nodes in the graph.
    steps : int, optional
        Number of logarithmic time samples (default 600).
    t1 : int, optional
        Lower exponent for the log-spaced time grid (default -2).
    t2 : int, optional
        Upper exponent for the log-spaced time grid (default 5).
    typf : type, optional
        Floating-point type used for computations (default np.float64).
    entropy_norm : str, optional
        Entropy normalization mode: "standard" or "complement".
    specific_heat_scale : str, optional
        Scaling for the entropy derivative: "logN" or "none".
    lanczos_steps : int, optional
        Number of Lanczos iterations (default 50).
    num_samples : int, optional
        Number of random probe vectors (default 30).
    seed : int, optional
        Random seed for reproducibility.
    tau_cutoff : float, optional
        Manual cutoff for numerical stability. If None, estimates from spectral gap.
    auto_cutoff_factor : float, default 10.0
        If tau_cutoff is None, sets cutoff to auto_cutoff_factor / lambda_min.
        Larger values = more aggressive cutoff. Typical: 10-100.

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray, NDArray]
        Normalized entropy profile, entropy derivative, variance profile,
        and the sampled time grid.
    """
    if laplacian.shape[0] != laplacian.shape[1]:
        raise ValueError("laplacian must be a square matrix.")
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")
    if num_samples < 1:
        raise ValueError("num_samples must be at least 1.")
    if lanczos_steps < 2:
        raise ValueError("lanczos_steps must be at least 2.")

    laplacian = (
        laplacian.astype(typf)
        if hasattr(laplacian, "astype")
        else np.asarray(laplacian, dtype=typf)
    )

    time_grid = np.logspace(t1, t2, steps, dtype=typf)

    # Determine tau cutoff (either manual or auto-estimated)
    if tau_cutoff is None:
        # Estimate spectral gap to set intelligent cutoff
        try:
            import scipy.sparse
            from scipy.sparse.linalg import eigsh

            # Ensure sparse format
            if not scipy.sparse.issparse(laplacian):
                L_sparse = scipy.sparse.csr_matrix(laplacian)
            else:
                L_sparse = laplacian

            # Compute smallest 10 eigenvalues (fast even for large matrices)
            n_eigen = min(10, laplacian.shape[0] - 2)
            eigenvals_small = eigsh(
                L_sparse, k=n_eigen, which="SM", return_eigenvectors=False
            )

            # Find spectral gap (smallest non-zero eigenvalue)
            threshold = 100 * np.finfo(typf).eps * np.max(np.abs(eigenvals_small))
            lambda_min = np.min(eigenvals_small[eigenvals_small > threshold])

            # Set cutoff: diffusion time ~ 1/lambda_min
            tau_cutoff = auto_cutoff_factor / lambda_min

        except Exception:
            # Fallback: heuristic based on graph size
            tau_cutoff = typf(num_nodes**1.5)

    # Identify tau values within stable range
    stable_mask = time_grid <= tau_cutoff

    # Only compute SLQ for stable tau values
    if np.any(stable_mask):
        time_grid_stable = time_grid[stable_mask]
        trace_exp, trace_l_exp, trace_l2_exp = _slq_trace_estimates(
            laplacian=laplacian,
            time_grid=time_grid_stable,
            num_samples=num_samples,
            lanczos_steps=lanczos_steps,
            seed=seed,
            typf=typf,
        )
    else:
        trace_exp = np.array([], dtype=typf)
        trace_l_exp = np.array([], dtype=typf)
        trace_l2_exp = np.array([], dtype=typf)

    log_N = np.log(typf(num_nodes)) if num_nodes > 0 else typf(1)
    entropy_profile = np.zeros_like(time_grid, dtype=typf)
    variance_profile = np.zeros_like(time_grid, dtype=typf)

    # Process stable tau values
    if np.any(stable_mask):
        valid_stable = trace_exp > 0
        avg_stable = np.zeros_like(trace_exp, dtype=typf)
        avg2_stable = np.zeros_like(trace_exp, dtype=typf)
        avg_stable[valid_stable] = trace_l_exp[valid_stable] / trace_exp[valid_stable]
        avg2_stable[valid_stable] = trace_l2_exp[valid_stable] / trace_exp[valid_stable]

        with np.errstate(divide="ignore", invalid="ignore"):
            entropy_profile[stable_mask] = np.where(
                valid_stable,
                (time_grid_stable * avg_stable + np.log(trace_exp)) / log_N,
                0.0,
            )
            variance_profile[stable_mask] = np.where(
                valid_stable,
                avg2_stable - avg_stable**2,
                0.0,
            )

    # For large tau (> cutoff), use analytical limiting values
    unstable_mask = ~stable_mask
    if np.any(unstable_mask):
        entropy_profile[unstable_mask] = log_N / log_N  # = 1.0
        variance_profile[unstable_mask] = 0.0

    normalized_entropy = _normalize_entropy_profile(entropy_profile, entropy_norm)

    # Compute specific heat
    entropy_derivative = np.gradient(normalized_entropy, np.log(time_grid))

    if np.any(unstable_mask):
        entropy_derivative[unstable_mask] = 0.0

    scale = specific_heat_scale.lower()
    if scale == "logn":
        entropy_derivative = log_N * entropy_derivative
    elif scale != "none":
        raise ValueError("specific_heat_scale must be 'logN' or 'none'.")

    return normalized_entropy, entropy_derivative, variance_profile, time_grid
