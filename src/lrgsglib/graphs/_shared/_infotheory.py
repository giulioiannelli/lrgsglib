"""Engine-agnostic information theory mixin functions.

These functions are imported into both SignedGraphNX and SignedGraphGT class
bodies. They operate on ``self.eigv``, ``self.N``, ``self.slp`` — attributes
present on both engines.
"""

import numpy as np
from typing import Optional
from numpy.typing import NDArray

from ...utils.basic import dtype_numerical_precision
from ...utils.lrg.infocomm import (
    entropy,
    compute_entropy_observables_slq,
    compute_renyi_observables_from_eigenvalues,
)
from ._backend import BackendManager, Backend


def _ensure_spectrum(self, typf=np.float64, backend="numpy"):
    """Ensure ``self.eigv`` is populated; compute spectrum if missing."""
    has_spectrum = hasattr(self, "eigv") and self.eigv is not None
    if has_spectrum:
        eigv_count = len(self.eigv)
        is_sparse = getattr(self, "_sparse_spectrum", False)
        expected_min = self.N - 2 if is_sparse else self.N
        if eigv_count >= expected_min:
            return  # already have enough eigenvalues

    # Both NX and GT expose compute_laplacian_spectrum or get_laplacian_spectrum
    if hasattr(self, "compute_laplacian_spectrum"):
        self.compute_laplacian_spectrum(typf=typf, backend=backend)
    elif hasattr(self, "get_laplacian_spectrum"):
        self.eigv = self.get_laplacian_spectrum(backend=backend)
    else:
        raise AttributeError(
            "Graph has no compute_laplacian_spectrum or get_laplacian_spectrum method."
        )


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
    entropy_mode: str = "exact",
    lanczos_steps: int = 50,
    num_samples: int = 30,
    slq_seed: Optional[int] = None,
    tau_cutoff: Optional[float] = None,
    auto_cutoff_factor: float = 10.0,
) -> None:
    """Compute the Shannon entropy profile of the signed Laplacian spectrum.

    Builds the normalized entropy H(tau) on a logarithmic time grid and its
    derivative ("specific heat"). In ``entropy_mode='exact'`` uses eigenvalues;
    in ``entropy_mode='slq'`` uses stochastic Lanczos quadrature.

    Results are stored as attributes: ``entropy``, ``specific_heat``,
    ``variance_profile``, ``tauscale``, ``entropy_params``.
    """
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    if backend is None:
        backend = getattr(self, "_backend_name", Backend.NUMPY.value)

    backend_obj = BackendManager.get_backend(backend, fallback=True)
    resolved_backend = backend_obj.name

    entropy_mode = entropy_mode.lower()
    if entropy_mode not in {"exact", "slq"}:
        raise ValueError("entropy_mode must be 'exact' or 'slq'.")

    if entropy_mode == "exact":
        _ensure_spectrum(self, typf=typf, backend=resolved_backend)

        eigv_count = len(self.eigv)
        is_sparse = getattr(self, "_sparse_spectrum", False)
        if is_sparse and eigv_count < self.N:
            import warnings
            warnings.warn(
                f"Using sparse spectrum with {eigv_count} eigenvalues "
                f"(missing {self.N - eigv_count}). Entropy may be less accurate.",
                UserWarning,
            )

    threshold = w_thresh if w_thresh is not None else dtype_numerical_precision(typf)
    backend_requested = backend

    if entropy_mode == "exact":
        eigenvalues = np.asarray(self.eigv, dtype=typf)
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
    else:
        if backend_requested == "cupy":
            import warnings
            warnings.warn(
                "SLQ entropy currently runs on CPU; ignoring backend='cupy'.",
                RuntimeWarning,
                stacklevel=2,
            )
        laplacian = self.slp.astype(typf)
        normalized_entropy, entropy_derivative, variance_profile, time_grid = (
            compute_entropy_observables_slq(
                laplacian=laplacian,
                num_nodes=self.N,
                steps=steps,
                t1=t1,
                t2=t2,
                typf=typf,
                entropy_norm=entropy_norm,
                specific_heat_scale=specific_heat_scale,
                lanczos_steps=lanczos_steps,
                num_samples=num_samples,
                seed=slq_seed,
                tau_cutoff=tau_cutoff,
                auto_cutoff_factor=auto_cutoff_factor,
            )
        )

    self.entropy = normalized_entropy
    self.specific_heat = entropy_derivative
    self.entropy_derivative = entropy_derivative  # deprecated alias
    self.variance_profile = variance_profile
    self.tauscale = time_grid
    self._last_entropy_backend = resolved_backend if entropy_mode == "exact" else "slq"
    self._last_entropy_request = backend_requested
    self.entropy_params = {
        "steps": steps,
        "t1": t1,
        "t2": t2,
        "w_thresh": threshold,
        "typf": typf,
        "backend": backend,
        "backend_requested": backend_requested,
        "backend_used": resolved_backend if entropy_mode == "exact" else "slq",
        "transpose": transpose,
        "entropy_norm": entropy_norm,
        "specific_heat_scale": specific_heat_scale,
        "entropy_mode": entropy_mode,
        "lanczos_steps": lanczos_steps,
        "num_samples": num_samples,
        "slq_seed": slq_seed,
    }


def get_entropy(self) -> NDArray:
    """Get cached entropy; compute with defaults if missing."""
    if not hasattr(self, "entropy") or self.entropy is None:
        compute_signed_laplacian_entropy(self)
    return self.entropy


def get_specific_heat(self) -> NDArray:
    """Get specific heat (entropy derivative)."""
    if not hasattr(self, "specific_heat") or self.specific_heat is None:
        compute_signed_laplacian_entropy(self)
    return self.specific_heat


def get_entropy_derivative(self) -> NDArray:
    """Deprecated: use ``get_specific_heat`` instead."""
    import warnings
    warnings.warn(
        "get_entropy_derivative() is deprecated, use get_specific_heat() instead",
        DeprecationWarning,
        stacklevel=2,
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
    """Compute Renyi/Tsallis entropy and generalized specific heat for ``q``.

    Results stored in ``self.renyi_results[q]`` and ``self.renyi_last``.
    """
    if q <= 0:
        raise ValueError("q must be strictly positive.")
    if steps < 2:
        raise ValueError("steps must be at least 2 to build the entropy profile.")
    if not (0 < tail_fraction <= 1):
        raise ValueError("tail_fraction must be in the interval (0, 1].")
    _ = transpose

    if backend is None:
        backend = getattr(self, "_backend_name", Backend.NUMPY.value)

    backend_obj = BackendManager.get_backend(backend, fallback=True)
    resolved_backend = backend_obj.name

    _ensure_spectrum(self, typf=typf, backend=resolved_backend)

    eigv_count = len(self.eigv)
    is_sparse = getattr(self, "_sparse_spectrum", False)
    if is_sparse and eigv_count < self.N:
        import warnings
        warnings.warn(
            f"Using sparse spectrum with {eigv_count} eigenvalues "
            f"(missing {self.N - eigv_count}). Entropy may be less accurate.",
            UserWarning,
        )

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
    """Retrieve cached Renyi profile for ``q``; compute with defaults if missing."""
    if (
        not hasattr(self, "renyi_results")
        or self.renyi_results is None
        or q not in self.renyi_results
    ):
        compute_renyi_entropy_profile(self, q=q)
    return self.renyi_results[q]
