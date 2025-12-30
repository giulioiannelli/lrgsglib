import numpy as np
from typing import Optional, TYPE_CHECKING
from numpy.typing import NDArray

from ...utils.basic import dtype_numerical_precision
from ...utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

# reuse spectral building blocks
from ._spectral import compute_laplacian_spectrum
from ._backend import BackendManager, Backend


if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def compute_signed_laplacian_entropy(
    self,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    w_thresh: float = None,
    typf: type = np.float64,
    backend: Optional[str] = None,
    transpose: Optional[bool] = None,
) -> None:
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    spectrum_missing = (
        not hasattr(self, "eigv") or self.eigv is None or len(self.eigv) != self.N
    )

    # Resolve backend using BackendManager with proper defaults and fallback
    if backend is None:
        resolved_backend = Backend.NUMPY.value
    else:
        backend_obj = BackendManager.get_backend(backend, fallback=True)
        resolved_backend = backend_obj.name

    if spectrum_missing:
        compute_laplacian_spectrum(self, typf=typf, backend=resolved_backend)

    eigenvalues = np.asarray(self.eigv, dtype=typf)
    threshold = w_thresh if w_thresh is not None else dtype_numerical_precision(typf)

    normalized_entropy, entropy_derivative, variance_profile, time_grid = (
        compute_entropy_observables_from_eigenvalues(
            eigenvalues=eigenvalues,
            num_nodes=self.N,
            steps=steps,
            t1=t1,
            t2=t2,
            typf=typf,
            threshold=threshold,
            pad_last=False,
        )
    )

    self.entropy = normalized_entropy
    self.specific_heat = entropy_derivative  # PRIMARY NAME
    self.entropy_derivative = entropy_derivative  # DEPRECATED ALIAS
    self.variance_profile = variance_profile
    self.tauscale = time_grid
    self.entropy_params = {
        "steps": steps,
        "t1": t1,
        "t2": t2,
        "w_thresh": threshold,
        "typf": typf,
        "backend": backend,
        "transpose": transpose,
    }


def get_entropy(self: "SignedGraph") -> NDArray:
    if not hasattr(self, "entropy") or self.entropy is None:
        compute_signed_laplacian_entropy(self)
    return self.entropy


def get_specific_heat(self: "SignedGraph") -> NDArray:
    """Get specific heat (entropy derivative)."""
    if not hasattr(self, "specific_heat") or self.specific_heat is None:
        compute_signed_laplacian_entropy(self)
    return self.specific_heat


def get_entropy_derivative(self: "SignedGraph") -> NDArray:
    """
    Get entropy derivative.

    .. deprecated::
        Use :func:`get_specific_heat` instead. This method is kept for
        backward compatibility.
    """
    import warnings
    warnings.warn(
        "get_entropy_derivative() is deprecated, use get_specific_heat() instead",
        DeprecationWarning,
        stacklevel=2
    )
    return self.get_specific_heat()


def compute_renyi_entropy_profile(
    self,
    q: float,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    w_thresh: float | None = None,
    typf: type = np.float64,
    backend: Optional[str] = None,
    transpose: Optional[bool] = None,
    tail_fraction: float = 0.1,
    drop_invalid: bool = True,
) -> None:
    """
    Compute Rényi/Tsallis entropy and generalized specific heat for a given ``q``.

    Laplacian eigenvalues play the role of energy levels. Entropies are
    normalized by ``log(N)`` to mirror the standard Shannon profile. For
    ``q→1`` this reduces to the normalized Shannon entropy.

    Note: eigenvectors are not used; ``transpose`` is kept for API compatibility
    and ignored.
    """
    if q <= 0:
        raise ValueError("q must be strictly positive.")
    if steps < 2:
        raise ValueError("steps must be at least 2 to build the entropy profile.")
    if not (0 < tail_fraction <= 1):
        raise ValueError("tail_fraction must be in the interval (0, 1].")
    _ = transpose  # kept for signature compatibility; eigenvectors not used

    spectrum_missing = (
        not hasattr(self, "eigv") or self.eigv is None or len(self.eigv) != self.N
    )

    # Resolve backend using BackendManager with proper defaults and fallback
    if backend is None:
        resolved_backend = Backend.NUMPY.value
    else:
        backend_obj = BackendManager.get_backend(backend, fallback=True)
        resolved_backend = backend_obj.name

    if spectrum_missing:
        compute_laplacian_spectrum(self, typf=typf, backend=resolved_backend)

    eigvals = np.asarray(self.eigv, dtype=typf)
    eps = w_thresh if w_thresh is not None else dtype_numerical_precision(typf)
    eigvals = np.where(np.abs(eigvals) > eps, eigvals, typf(0))
    log_N = np.log(typf(self.N)) if self.N > 1 else typf(1)

    tau_grid = np.logspace(t1, t2, steps, dtype=typf)
    renyi_entropy = np.full_like(tau_grid, np.nan, dtype=typf)
    tsallis_entropy = np.full_like(tau_grid, np.nan, dtype=typf)
    log_tau = np.log(tau_grid)

    invalid_counts = 0
    for idx, tau in enumerate(tau_grid):
        if abs(q - 1.0) < 1e-8:
            weights = np.exp(-tau * eigvals)
            Z = np.nansum(weights, dtype=typf)
            if Z == 0:
                continue
            p = weights / Z
            with np.errstate(divide="ignore", invalid="ignore"):
                shannon_raw = -np.nansum(p * np.log(p), dtype=typf)
            shannon_norm = shannon_raw / log_N
            renyi_entropy[idx] = shannon_norm
            tsallis_entropy[idx] = shannon_norm
            continue

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
        renyi_entropy[idx] = renyi_raw / log_N
        tsallis_entropy[idx] = tsallis_raw / log_N

    with np.errstate(invalid="ignore"):
        specific_heat = -np.gradient(renyi_entropy, log_tau)

    tail_len = max(1, int(round(tail_fraction * steps)))
    tail_slice = specific_heat[-tail_len:]
    ds_estimate = (
        np.nan if tail_slice.size == 0 else typf(2.0) * np.nanmean(tail_slice, dtype=typf)
    )

    results = {
        "q": q,
        "tau": tau_grid,
        "renyi_entropy": renyi_entropy,
        "tsallis_entropy": tsallis_entropy,
        "specific_heat": specific_heat,
        "ds_estimate": ds_estimate,
        "tail_fraction": tail_fraction,
        "invalid_weight_counts": invalid_counts,
    }

    if not hasattr(self, "renyi_results") or self.renyi_results is None:
        self.renyi_results = {}
    self.renyi_results[q] = results
    self.renyi_last = results


def get_renyi_results(self, q: float) -> dict:
    """
    Retrieve cached Rényi profile for ``q``; compute with defaults if missing.
    """
    if not hasattr(self, "renyi_results") or self.renyi_results is None or q not in self.renyi_results:
        compute_renyi_entropy_profile(self, q=q)
    return self.renyi_results[q]
