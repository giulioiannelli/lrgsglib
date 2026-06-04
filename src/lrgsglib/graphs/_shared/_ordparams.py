"""Engine-agnostic order parameter mixin functions (spectral gap)."""

import numpy as np
from typing import Optional

from ...config.const import SG_LAPL_DEFAULT_TYPE


def compute_gap(
    self,
    backend: Optional[str] = None,
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
    laplacian_type: str = SG_LAPL_DEFAULT_TYPE,
) -> None:
    """Compute and cache the spectral gap (eigenvalue 1 - eigenvalue 0)."""
    compute_gap_between(
        self,
        low=0,
        high=1,
        backend=backend,
        typf=typf,
        transpose=transpose,
        flip_to_pos=flip_to_pos,
        laplacian_type=laplacian_type,
    )


def compute_gap_between(
    self,
    low: int = 0,
    high: int = 1,
    backend: Optional[str] = None,
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
    rescale_by_sqrt: bool = True,
    laplacian_type: str = SG_LAPL_DEFAULT_TYPE,
) -> None:
    """Compute gap between ``eigv[low]`` and ``eigv[high]`` of the requested
    Laplacian (``laplacian_type`` in {'signed','sym','rw'})."""
    need_compute = (
        not hasattr(self, "eigv")
        or self.eigv is None
        or self.eigv.size != self.N
        or getattr(self, "_spectrum_laplacian_type", None) != laplacian_type
    )
    if need_compute:
        kw = {
            "typf": typf,
            "transpose": transpose,
            "flip_to_pos": flip_to_pos,
            "laplacian_type": laplacian_type,
        }
        if backend is not None:
            kw["backend"] = backend
        # Filter kwargs to only pass what the method accepts
        import inspect
        if hasattr(self, "compute_laplacian_spectrum_weigV"):
            sig = inspect.signature(self.compute_laplacian_spectrum_weigV)
            valid_kw = {k: v for k, v in kw.items() if k in sig.parameters}
            self.compute_laplacian_spectrum_weigV(**valid_kw)
        else:
            raise AttributeError("No method to compute eigendecomposition.")

    eigv = self.eigv
    if eigv is None or eigv.size < 2:
        raise ValueError("Cannot compute gap: insufficient eigenvalues.")
    if not (isinstance(low, int) and isinstance(high, int)):
        raise TypeError("low and high must be integer indices")
    N = eigv.size
    if not (0 <= low < high < N):
        raise ValueError(
            f"Indices must satisfy 0 <= low < high < N ({N}); "
            f"got low={low}, high={high}"
        )

    diff = float(eigv[high] - eigv[low])
    largest = float(eigv[-1])
    gap = diff / largest if largest != 0 else diff
    if rescale_by_sqrt:
        gap *= np.sqrt(self.N)
    self.gap = gap
    self._gap_laplacian_type = laplacian_type


def get_gap(self, laplacian_type: str = SG_LAPL_DEFAULT_TYPE) -> float:
    """Return cached gap for ``laplacian_type``; compute it if missing/stale."""
    if (
        not hasattr(self, "gap")
        or self.gap is None
        or getattr(self, "_gap_laplacian_type", None) != laplacian_type
    ):
        compute_gap(self, laplacian_type=laplacian_type)
    return self.gap
