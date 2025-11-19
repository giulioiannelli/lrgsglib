import numpy as np

from typing import Optional, TYPE_CHECKING
from lrgsglib.config.const import SG_REPR


if TYPE_CHECKING:
    # avoid runtime circular import; used only for type checking / annotations
    from .SignedGraph import SignedGraph

# import spectral helpers directly so functions here can call them
# as plain functions with `self` as first argument (avoids relying on
# class-level bound imports on SignedGraph)
from ._spectral import compute_laplacian_spectrum_weigV


def compute_gap(
    self: "SignedGraph",
    backend: str = 'numpy',
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
    low: int = 0,
    high: int = 1,
) -> float:
    """Compute and cache the spectral gap of the signed Laplacian for `self`.

    This forwards to `compute_gap_between` which accepts arbitrary
    eigenvalue indices. By default it computes the low=0 / high=1 gap
    (same behaviour as before).
    """
    # forward to the more general implementation which handles arbitrary
    # eigenvalue indices; compute_gap_between will validate indices.
    gap = compute_gap_between(
        self,
        low=low,
        high=high,
        backend=backend,
        typf=typf,
        transpose=transpose,
        flip_to_pos=flip_to_pos,
    )

    # cache the default (low=0, high=1) gap under _gap for compatibility
    if low == 0 and high == 1:
        self._gap = gap
    return gap


def compute_gap_between(
    self: "SignedGraph",
    low: int = 0,
    high: int = 1,
    backend: str = 'numpy',
    typf: type = np.float64,
    transpose: bool = True,
    flip_to_pos: bool = True,
) -> float:
    """Compute gap between eigenvalues `eigv[low]` and `eigv[high]`.

    Parameters
    ----------
    low, high : int
        Indices of the eigenvalues to compute the gap between. Must satisfy
        `0 <= low < high < N` where `N` is the number of eigenvalues.
    Other args : forwarded to spectrum computation if needed.

    Returns
    -------
    float
        The computed (possibly normalized) gap.
    """
    # Ensure eigenvalues are present and complete
    if not hasattr(self, 'eigv') or self.eigv is None or self.eigv.size != self.N:
        compute_laplacian_spectrum_weigV(
            self,
            typf=typf,
            transpose=transpose,
            backend=backend,
            flip_to_pos=flip_to_pos,
        )

    eigv = self.eigv
    if eigv is None or eigv.size < 2:
        raise ValueError('Cannot compute gap: insufficient eigenvalues.')

    # validate indices
    if not (isinstance(low, int) and isinstance(high, int)):
        raise TypeError('low and high must be integer indices')
    N = eigv.size
    if not (0 <= low < high < N):
        raise ValueError(f'Indices must satisfy 0 <= low < high < N ({N}); got low={low}, high={high}')

    diff = float(eigv[high] - eigv[low])
    largest = float(eigv[-1])
    if largest == 0:
        gap = diff
    else:
        gap = diff / largest

    return float(gap)


def get_gap(self: "SignedGraph") -> float:
    """Return cached gap; compute it if missing."""
    if hasattr(self, '_gap') and self._gap is not None:
        return self._gap
    return compute_gap(self)


def compute_pinf(
    self: "SignedGraph",
    which: int = 0,
    val: Optional[object] = 1,
    on_g: str = SG_REPR,
) -> None:
    """Compute the infinite-cluster probability (Pinf) for `self`.

    Mirrors the former method implementation attached to `SignedGraph`.
    """
    if val is None:
        val = 1

    # Delegate to the SignedGraph helper for cluster sizes
    clustd = np.array(self.get_eigV_cluster_sizes(which, True, val, on_g))
    self.Pinf = clustd[0] / self.N

    denominator = np.sum(clustd) - clustd[0]
    if denominator > 0:
        self.Pinf_var = np.sum(clustd @ clustd - clustd[0] ** 2) / denominator
    else:
        self.Pinf_var = 0.0

    if hasattr(self, "Pinf_dict"):
        self.Pinf_dict[which] = (self.Pinf, self.Pinf_var)
    else:
        self.Pinf_dict = {which: (self.Pinf, self.Pinf_var)}
