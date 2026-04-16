"""Classical entropy computation from eigenvalues.

This module provides the standard eigenvalue-based entropy computation
using the von Neumann entropy formula.
"""

from typing import Callable, Optional, Tuple

import numpy as np
from networkx import Graph
from numpy.typing import NDArray

from ...basic import dtype_numerical_precision
from ..spectral import get_graph_lspectrum
from ._common import _normalize_entropy_profile

__all__ = [
    "entropy",
    "compute_entropy_observables_from_eigenvalues",
    "specific_heat_tau_window",
]


def specific_heat_tau_window(
    tau_grid: NDArray,
    specific_heat: NDArray,
    tau_mark: Optional[float] = None,
    rel_threshold: float = 0.01,
    log_pad: float = 0.5,
) -> Tuple[float, float]:
    """Return a (``tau_lo``, ``tau_hi``) window covering the meaningful range of
    the LRG specific heat ``C(tau)``.

    The window is the sub-range where ``C(tau) > rel_threshold * C(tau*)``
    (with ``tau* = argmax C``), padded by ``log_pad`` decades on each side
    and guaranteed to contain ``tau_mark`` when supplied. Intended for
    trimming ``semilogx`` plots of LRG observables so the interesting
    crossover region is shown without wasted empty space on either side.

    Parameters
    ----------
    tau_grid : NDArray
        Strictly monotonically increasing array of RG scales on which
        ``specific_heat`` was evaluated.
    specific_heat : NDArray
        The specific-heat profile ``C(tau)``, same length as ``tau_grid``.
    tau_mark : float, optional
        A reference scale (e.g. ``1/lambda_max``) that must stay inside the
        returned window. If ``None``, no such constraint is applied.
    rel_threshold : float, optional
        Fraction of the peak height used as the significance cut-off
        (default ``0.01`` = 1 % of the peak).
    log_pad : float, optional
        Decades of padding added on each side in ``log10(tau)``
        (default ``0.5``).

    Returns
    -------
    Tuple[float, float]
        ``(tau_lo, tau_hi)`` clamped to ``[tau_grid[0], tau_grid[-1]]``.
    """
    tau_grid = np.asarray(tau_grid, dtype=float)
    C = np.asarray(specific_heat, dtype=float)
    if tau_grid.shape != C.shape:
        raise ValueError("tau_grid and specific_heat must have the same shape.")
    if tau_grid.size == 0:
        raise ValueError("tau_grid must be non-empty.")

    peak = int(np.argmax(C))
    mask = C > rel_threshold * C[peak]
    if not mask.any():
        return float(tau_grid[0]), float(tau_grid[-1])

    idx = np.where(mask)[0]
    lo = float(tau_grid[idx[0]]) * 10.0 ** (-log_pad)
    hi = float(tau_grid[idx[-1]]) * 10.0 ** (+log_pad)
    if tau_mark is not None:
        lo = min(lo, float(tau_mark) * 10.0 ** (-log_pad))
        hi = max(hi, float(tau_mark) * 10.0 ** (+log_pad))
    return max(lo, float(tau_grid[0])), min(hi, float(tau_grid[-1]))


def entropy(
    G: Optional[Graph] = None,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    wTresh: float = dtype_numerical_precision(np.float64),
    func_lspectrum: Optional[Callable] = get_graph_lspectrum,
    eigenvalues: Optional[NDArray] = None,
    num_nodes: Optional[int] = None,
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
    typf: type = np.float64,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Compute the entropy, its derivative, variance, and time steps for a graph.

    Parameters
    ----------
    G : networkx.Graph, optional
        The input graph. Required if ``eigenvalues`` is not provided.
    steps : int, optional
        Number of time steps for the computation (default is 600).
    t1 : int, optional
        Logarithmic start time (default is -2).
    t2 : int, optional
        Logarithmic end time (default is 5).
    wTresh : float, optional
        Threshold for filtering eigenvalues (default is machine precision for
        float64). Use None to infer from ``typf``.
    func_lspectrum : Callable, optional
        Function to compute the Laplacian spectrum of the graph (default is
        `get_graph_lspectrum`). Required if ``eigenvalues`` is not provided.
    eigenvalues : NDArray, optional
        Precomputed Laplacian spectrum. If provided, ``G`` is optional.
    num_nodes : int, optional
        Number of nodes for entropy normalization. Defaults to ``len(eigenvalues)``
        when ``eigenvalues`` is provided.
    entropy_norm : str, optional
        Entropy normalization mode: "standard" or "complement".
    specific_heat_scale : str, optional
        Scaling for the entropy derivative: "logN" or "none".
    typf : type, optional
        Floating-point type used for computations (default np.float64).

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray, NDArray]
        - Normalized entropy values.
        - Derivative of entropy with respect to time.
        - Variance of the Laplacian spectrum.
        - Time steps used for the computation.
    """
    if wTresh is None:
        wTresh = dtype_numerical_precision(typf)

    if eigenvalues is None:
        if G is None:
            raise ValueError("Provide G or eigenvalues to compute entropy.")
        if func_lspectrum is None:
            raise ValueError("func_lspectrum must be provided when G is used.")
        N = G.number_of_nodes()
        _, w = func_lspectrum(G)
    else:
        w = eigenvalues
        N = num_nodes if num_nodes is not None else int(np.size(w))

    normalized_entropy, entropy_derivative, variance_profile, time_grid = (
        compute_entropy_observables_from_eigenvalues(
            eigenvalues=w,
            num_nodes=N,
            steps=steps,
            t1=t1,
            t2=t2,
            typf=typf,
            threshold=wTresh,
            entropy_norm=entropy_norm,
            specific_heat_scale=specific_heat_scale,
        )
    )

    return normalized_entropy, entropy_derivative, variance_profile, time_grid


def compute_entropy_observables_from_eigenvalues(
    eigenvalues: NDArray,
    num_nodes: int,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    typf: type = np.float64,
    threshold: Optional[float] = None,
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Compute entropy-related observables from a Laplacian spectrum.

    Parameters
    ----------
    eigenvalues : NDArray
        Spectrum of the Laplacian (or signed Laplacian) matrix.
    num_nodes : int
        Number of nodes in the graph used to normalize the entropy.
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
    entropy_norm : str, optional
        Entropy normalization mode: "standard" or "complement".
    specific_heat_scale : str, optional
        Scaling for the entropy derivative: "logN" or "none".

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray, NDArray]
        Normalized entropy profile, entropy derivative, variance profile,
        and the sampled time grid.
    """
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    eigvals = np.asarray(eigenvalues, dtype=typf)
    eps = threshold if threshold is not None else dtype_numerical_precision(typf)
    eigvals = np.where(np.abs(eigvals) > eps, eigvals, typf(0))

    time_grid = np.logspace(t1, t2, steps, dtype=typf)
    entropy_profile = np.zeros_like(time_grid, dtype=typf)
    variance_profile = np.zeros_like(time_grid, dtype=typf)

    log_N = np.log(typf(num_nodes)) if num_nodes > 0 else typf(1)

    for idx, tau in enumerate(time_grid):
        rhoTr = np.exp(-tau * eigvals)
        trace_rho = np.nansum(rhoTr, dtype=typf)
        if trace_rho == 0:
            rho = np.zeros_like(rhoTr)
        else:
            rho = rhoTr / trace_rho

        with np.errstate(divide="ignore", invalid="ignore"):
            entropy_profile[idx] = -np.nansum(rho * np.log(rho), dtype=typf) / log_N

        if trace_rho:
            avgrho = np.nansum(eigvals * rhoTr, dtype=typf) / trace_rho
            av2rho = np.nansum((eigvals**2) * rhoTr, dtype=typf) / trace_rho
        else:
            avgrho = typf(0)
            av2rho = typf(0)
        variance_profile[idx] = av2rho - avgrho**2

    normalized_entropy = _normalize_entropy_profile(entropy_profile, entropy_norm)

    # Specific heat is always computed from the complement form (1 - S/logN)
    # so that C(τ) = d(1-S/logN)/d(logτ) stays positive at diffusion scales,
    # regardless of which normalization the caller chose for the output entropy.
    complement = _normalize_entropy_profile(entropy_profile, "complement")
    entropy_derivative = np.gradient(complement, np.log(time_grid))
    scale = specific_heat_scale.lower()
    if scale == "logn":
        entropy_derivative = log_N * entropy_derivative
    elif scale != "none":
        raise ValueError("specific_heat_scale must be 'logN' or 'none'.")
    return normalized_entropy, entropy_derivative, variance_profile, time_grid
