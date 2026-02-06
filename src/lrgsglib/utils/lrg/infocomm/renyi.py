"""Rényi and Tsallis entropy computation.

This module provides generalized entropy computation using Rényi and
Tsallis entropy formulas with parameterized q values.
"""

from typing import Any, Dict, Optional

import numpy as np
from numpy.typing import NDArray

from ...basic import dtype_numerical_precision
from ._common import _normalize_entropy_profile

__all__ = [
    "compute_renyi_observables_from_eigenvalues",
]


def compute_renyi_observables_from_eigenvalues(
    eigenvalues: NDArray,
    num_nodes: int,
    q: float,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    typf: type = np.float64,
    threshold: Optional[float] = None,
    tail_fraction: float = 0.1,
    drop_invalid: bool = True,
    entropy_norm: Optional[str] = None,
    legacy_normalization: bool = False,
    specific_heat_scale: str = "none",
) -> Dict[str, Any]:
    """
    Compute Renyi/Tsallis entropy profiles and generalized specific heat.

    Parameters
    ----------
    eigenvalues : NDArray
        Spectrum of the Laplacian (or signed Laplacian) matrix.
    num_nodes : int
        Number of nodes in the graph used to normalize the entropy.
    q : float
        Order of the Renyi/Tsallis entropy (must be > 0).
    steps : int, optional
        Number of logarithmic time samples (default 600).
    t1 : int, optional
        Lower exponent for the log-spaced time grid (default -2).
    t2 : int, optional
        Upper exponent for the log-spaced time grid (default 5).
    typf : type, optional
        Floating-point type used for computations (default np.float64).
    threshold : float, optional
        Magnitude threshold below which eigenvalues are treated as zero.
        Defaults to the machine precision of `typf`.
    tail_fraction : float, optional
        Fraction of the tail used for ds estimate (default 0.1).
    drop_invalid : bool, optional
        If True, discard invalid weights for q != 1 (default True).
    entropy_norm : str, optional
        Entropy normalization mode: "standard" or "complement". If None,
        defaults to "standard" in legacy mode and "complement" otherwise.
    legacy_normalization : bool, optional
        Compatibility mode for older results. When True and ``entropy_norm`` is
        not provided, the entropy profile is returned in "standard" form while
        the specific heat is computed from the complement.
    specific_heat_scale : str, optional
        Scaling for the entropy derivative: "logN" or "none".

    Returns
    -------
    dict
        Dictionary with keys:
        - "q": q value
        - "tau": time grid
        - "renyi_entropy": normalized Renyi entropy profile
        - "tsallis_entropy": normalized Tsallis entropy profile
        - "specific_heat": generalized specific heat profile
        - "ds_estimate": tail estimate of spectral dimension
        - "tail_fraction": tail fraction used
        - "invalid_weight_counts": count of invalid weights discarded
        - "entropy_norm": normalization mode used
        - "legacy_normalization": whether legacy normalization was used
        - "specific_heat_scale": scaling mode for the entropy derivative
    """
    if q <= 0:
        raise ValueError("q must be strictly positive.")
    if steps < 2:
        raise ValueError("steps must be at least 2 to build the entropy profile.")
    if not (0 < tail_fraction <= 1):
        raise ValueError("tail_fraction must be in the interval (0, 1].")

    eigvals = np.asarray(eigenvalues, dtype=typf)
    eps = threshold if threshold is not None else dtype_numerical_precision(typf)
    eigvals = np.where(np.abs(eigvals) > eps, eigvals, typf(0))
    log_N = np.log(typf(num_nodes)) if num_nodes > 1 else typf(1)

    if entropy_norm is None:
        entropy_norm = "standard" if legacy_normalization else "complement"

    tau_grid = np.logspace(t1, t2, steps, dtype=typf)
    renyi_entropy_standard = np.full_like(tau_grid, np.nan, dtype=typf)
    tsallis_entropy_standard = np.full_like(tau_grid, np.nan, dtype=typf)
    log_tau = np.log(tau_grid)

    invalid_counts = 0
    for idx, tau in enumerate(tau_grid):
        if abs(q - 1.0) < 1e-8:
            # Shannon entropy limit (q -> 1)
            weights = np.exp(-tau * eigvals)
            Z = np.nansum(weights, dtype=typf)
            if Z == 0:
                continue
            p = weights / Z
            with np.errstate(divide="ignore", invalid="ignore"):
                shannon_raw = -np.nansum(p * np.log(p), dtype=typf)
            shannon_norm = shannon_raw / log_N
            renyi_entropy_standard[idx] = shannon_norm
            tsallis_entropy_standard[idx] = shannon_norm
            continue

        # General q case
        base = 1.0 + (1.0 - q) * tau * eigvals
        if drop_invalid:
            valid = base > 0
            invalid_counts += np.count_nonzero(~valid)
            base = np.where(valid, base, 0.0)
        exponent = typf(1.0 / (q - 1.0))
        with np.errstate(over="ignore", invalid="ignore"):
            eps_i = np.power(base, exponent, dtype=typf)
        Z = np.nansum(eps_i, dtype=typf)
        if Z <= 0:
            continue
        p = eps_i / Z
        mq = np.nansum(np.power(p, q, dtype=typf), dtype=typf)
        if mq <= 0:
            continue
        renyi_raw = (1.0 / (1.0 - q)) * np.log(mq)
        tsallis_raw = (1.0 / (q - 1.0)) * (1.0 - mq)
        renyi_entropy_standard[idx] = renyi_raw / log_N
        tsallis_entropy_standard[idx] = tsallis_raw / log_N

    renyi_entropy = _normalize_entropy_profile(
        renyi_entropy_standard,
        entropy_norm,
    )
    tsallis_entropy = _normalize_entropy_profile(
        tsallis_entropy_standard,
        entropy_norm,
    )

    with np.errstate(invalid="ignore"):
        heat_mode = entropy_norm
        if legacy_normalization and entropy_norm.lower() == "standard":
            heat_mode = "complement"
        heat_source = _normalize_entropy_profile(
            renyi_entropy_standard,
            heat_mode,
        )
        specific_heat = np.gradient(heat_source, log_tau)
        scale = specific_heat_scale.lower()
        if scale == "logn":
            specific_heat = log_N * specific_heat
        elif scale != "none":
            raise ValueError("specific_heat_scale must be 'logN' or 'none'.")

    tail_len = max(1, int(round(tail_fraction * steps)))
    tail_slice = specific_heat[-tail_len:]
    ds_estimate = (
        np.nan
        if tail_slice.size == 0
        else typf(2.0) * np.nanmean(tail_slice, dtype=typf)
    )

    return {
        "q": q,
        "tau": tau_grid,
        "renyi_entropy": renyi_entropy,
        "tsallis_entropy": tsallis_entropy,
        "specific_heat": specific_heat,
        "ds_estimate": ds_estimate,
        "tail_fraction": tail_fraction,
        "invalid_weight_counts": invalid_counts,
        "entropy_norm": entropy_norm,
        "legacy_normalization": legacy_normalization,
        "specific_heat_scale": specific_heat_scale,
    }
