import numpy as np
from typing import Optional, TYPE_CHECKING
from numpy.typing import NDArray

from ...utils.basic import dtype_numerical_precision
from ...utils.lrg.infocomm import (
    entropy,
    compute_renyi_observables_from_eigenvalues,
)

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
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
) -> None:
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    spectrum_missing = (
        not hasattr(self, "eigv") or self.eigv is None or len(self.eigv) != self.N
    )

    # Resolve backend using instance default unless explicitly overridden.
    if backend is None:
        resolved_backend = getattr(self, "_backend_name", Backend.NUMPY.value)
    else:
        backend_obj = BackendManager.get_backend(backend, fallback=True)
        resolved_backend = backend_obj.name

    if spectrum_missing:
        compute_laplacian_spectrum(self, typf=typf, backend=resolved_backend)

    eigenvalues = np.asarray(self.eigv, dtype=typf)
    threshold = w_thresh if w_thresh is not None else dtype_numerical_precision(typf)

    normalized_entropy, entropy_derivative, variance_profile, time_grid = entropy(
        eigenvalues=eigenvalues,
        num_nodes=self.N,
        steps=steps,
        t1=t1,
        t2=t2,
        wTresh=threshold,
        entropy_norm=entropy_norm,
        specific_heat_scale=specific_heat_scale,
        typf=typf,
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
        "entropy_norm": entropy_norm,
        "specific_heat_scale": specific_heat_scale,
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
    entropy_norm: Optional[str] = None,
    legacy_normalization: bool = True,
    specific_heat_scale: str = "none",
) -> None:
    """
    Compute Rényi/Tsallis entropy and generalized specific heat for a given ``q``.

    Laplacian eigenvalues play the role of energy levels. Entropies are
    normalized by ``log(N)`` to mirror the standard Shannon profile. For
    ``q→1`` this reduces to the normalized Shannon entropy. Use
    ``entropy_norm="complement"`` to match the SignedGraph Shannon
    normalization; set ``legacy_normalization=True`` to preserve older outputs.
    Use ``specific_heat_scale="logN"`` to scale the derivative by ``log(N)``.

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

    # Resolve backend using instance default unless explicitly overridden.
    if backend is None:
        resolved_backend = getattr(self, "_backend_name", Backend.NUMPY.value)
    else:
        backend_obj = BackendManager.get_backend(backend, fallback=True)
        resolved_backend = backend_obj.name

    if spectrum_missing:
        compute_laplacian_spectrum(self, typf=typf, backend=resolved_backend)

    eigvals = np.asarray(self.eigv, dtype=typf)
    eps = w_thresh if w_thresh is not None else dtype_numerical_precision(typf)

    results = compute_renyi_observables_from_eigenvalues(
        eigenvalues=eigvals,
        num_nodes=self.N,
        q=q,
        steps=steps,
        t1=t1,
        t2=t2,
        typf=typf,
        threshold=eps,
        tail_fraction=tail_fraction,
        drop_invalid=drop_invalid,
        entropy_norm=entropy_norm,
        legacy_normalization=legacy_normalization,
        specific_heat_scale=specific_heat_scale,
    )

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
