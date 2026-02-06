import numpy as np
from typing import Optional, TYPE_CHECKING
from numpy.typing import NDArray

from ...utils.basic import dtype_numerical_precision
from ...utils.lrg.infocomm import (
    entropy,
    compute_entropy_observables_slq,
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
    entropy_mode: str = "exact",
    lanczos_steps: int = 50,
    num_samples: int = 30,
    slq_seed: Optional[int] = None,
    tau_cutoff: Optional[float] = None,
    auto_cutoff_factor: float = 10.0,
) -> None:
    """
    Compute the Shannon entropy profile of the signed Laplacian spectrum.

    Builds the normalized entropy H(τ) on a logarithmic time grid and its
    derivative ("specific heat"). In ``entropy_mode='exact'`` it uses the
    Laplacian eigenvalues (computing them if missing); in
    ``entropy_mode='slq'`` it estimates the observables via stochastic Lanczos
    quadrature (SLQ) directly from the Laplacian matrix.

    Parameters
    ----------
    steps : int, default 600
        Number of time samples between ``10**t1`` and ``10**t2``.
        Must be at least 1.
    t1 : int, default -2
        Base-10 exponent for the first time scale, i.e. ``τ_min = 10**t1``.
    t2 : int, default 5
        Base-10 exponent for the last time scale, i.e. ``τ_max = 10**t2``.
    w_thresh : float, optional
        Numerical threshold for near-zero eigenvalues. If ``None``, a value
        derived from the numerical precision of ``typf`` is used.
    typf : dtype, default numpy.float64
        Floating-point dtype used for computations.
    backend : {"numpy", "cupy"}, optional
        Backend for spectrum computation when ``entropy_mode='exact'``. If
        ``None``, uses the instance default backend.
    transpose : bool, optional
        Ignored. Kept for API compatibility with older versions.
    entropy_norm : {"complement", "standard"}, default "complement"
        Entropy normalization mode; forwarded to the internal ``entropy``
        routine.
    specific_heat_scale : {"logN", "none"}, default "logN"
        Scaling applied to the entropy derivative (specific heat);
        forwarded to the internal ``entropy`` routine.
    entropy_mode : {"exact", "slq"}, default "exact"
        Computation mode. ``'exact'`` uses eigenvalues (dense spectrum
        preferred; sparse with fewer than ``N`` eigenvalues is accepted but
        less accurate). ``'slq'`` uses stochastic Lanczos quadrature on CPU.
    lanczos_steps : int, default 50
        Number of Lanczos iterations used by SLQ when
        ``entropy_mode='slq'``.
    num_samples : int, default 30
        Number of random probe vectors used by SLQ.
    slq_seed : int, optional
        Random seed for SLQ probe generation.
    tau_cutoff : float, optional
        Manual cutoff for numerical stability in SLQ mode. If None (default),
        automatically estimates from spectral gap (recommended). For tau > tau_cutoff,
        uses analytical limiting values. Set to float('inf') to disable cutoff (risky).
    auto_cutoff_factor : float, default 10.0
        When tau_cutoff is None, sets cutoff to auto_cutoff_factor / lambda_min
        where lambda_min is the spectral gap. Larger = more conservative.
        Typical range: 10-100.

    Warns
    -----
    UserWarning
        When using a sparse spectrum with fewer than ``N`` eigenvalues in
        ``entropy_mode='exact'``; accuracy may be reduced.
    RuntimeWarning
        When ``backend='cupy'`` is requested together with
        ``entropy_mode='slq'``; SLQ currently runs on CPU.

    Returns
    -------
    None
        Results are stored as attributes on the instance:
        - ``entropy``: array of shape ``(steps,)``
        - ``specific_heat``: array of shape ``(steps,)``
        - ``entropy_derivative``: deprecated alias of ``specific_heat``
        - ``variance_profile``: array of shape ``(steps,)``
        - ``tauscale``: time grid of shape ``(steps,)``
        Parameters used are recorded in ``entropy_params``.

    Notes
    -----
    - If the spectrum is missing in ``entropy_mode='exact'``, it is computed
      via ``compute_laplacian_spectrum`` using the resolved backend.
    - A near-zero eigenvalue threshold is derived from ``typf`` via
      ``dtype_numerical_precision`` when ``w_thresh`` is not provided.
    """
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    # Resolve backend using instance default unless explicitly overridden
    backend_requested = backend
    if backend is None:
        backend = getattr(self, "_backend_name", Backend.NUMPY.value)

    backend_obj = BackendManager.get_backend(backend, fallback=True)
    resolved_backend = backend_obj.name

    entropy_mode = entropy_mode.lower()
    if entropy_mode not in {"exact", "slq"}:
        raise ValueError("entropy_mode must be 'exact' or 'slq'.")

    # Check if spectrum exists and has correct size
    # Accept N-2 eigenvalues if computed with sparse method
    has_spectrum = hasattr(self, "eigv") and self.eigv is not None
    if has_spectrum:
        eigv_count = len(self.eigv)
        is_sparse_spectrum = getattr(self, "_sparse_spectrum", False)

        # N-2 eigenvalues valid for sparse, N required for dense
        expected_min = self.N - 2 if is_sparse_spectrum else self.N
        spectrum_missing = eigv_count < expected_min
    else:
        spectrum_missing = True

    if entropy_mode == "exact" and spectrum_missing:
        compute_laplacian_spectrum(self, typf=typf, backend=resolved_backend)

    threshold = w_thresh if w_thresh is not None else dtype_numerical_precision(typf)

    if entropy_mode == "exact":
        # After computation, verify we have usable spectrum
        eigv_count = len(self.eigv)
        is_sparse_spectrum = getattr(self, "_sparse_spectrum", False)

        # Warn if using sparse spectrum for entropy (less accurate)
        if is_sparse_spectrum and eigv_count < self.N:
            import warnings
            warnings.warn(
                f"Using sparse spectrum with {eigv_count} eigenvalues (missing {self.N - eigv_count}). "
                f"Entropy calculation may be less accurate. Consider recomputing with dense method.",
                UserWarning
            )

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
    self.specific_heat = entropy_derivative  # PRIMARY NAME
    self.entropy_derivative = entropy_derivative  # DEPRECATED ALIAS
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

    # Resolve backend using instance default unless explicitly overridden
    if backend is None:
        backend = getattr(self, "_backend_name", Backend.NUMPY.value)

    backend_obj = BackendManager.get_backend(backend, fallback=True)
    resolved_backend = backend_obj.name

    # Check if spectrum exists and has correct size
    # Accept N-2 eigenvalues if computed with sparse method
    has_spectrum = hasattr(self, "eigv") and self.eigv is not None
    if has_spectrum:
        eigv_count = len(self.eigv)
        is_sparse_spectrum = getattr(self, "_sparse_spectrum", False)

        # N-2 eigenvalues valid for sparse, N required for dense
        expected_min = self.N - 2 if is_sparse_spectrum else self.N
        spectrum_missing = eigv_count < expected_min
    else:
        spectrum_missing = True

    if spectrum_missing:
        compute_laplacian_spectrum(self, typf=typf, backend=resolved_backend)

    # After computation, verify we have usable spectrum
    eigv_count = len(self.eigv)
    is_sparse_spectrum = getattr(self, "_sparse_spectrum", False)

    # Warn if using sparse spectrum for entropy (less accurate)
    if is_sparse_spectrum and eigv_count < self.N:
        import warnings
        warnings.warn(
            f"Using sparse spectrum with {eigv_count} eigenvalues (missing {self.N - eigv_count}). "
            f"Entropy calculation may be less accurate. Consider recomputing with dense method.",
            UserWarning
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
    """
    Retrieve cached Rényi profile for ``q``; compute with defaults if missing.
    """
    if not hasattr(self, "renyi_results") or self.renyi_results is None or q not in self.renyi_results:
        compute_renyi_entropy_profile(self, q=q)
    return self.renyi_results[q]
