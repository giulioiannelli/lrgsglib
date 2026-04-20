"""Offline overlap kernels for signed-walker visit-count fields.

All functions are pure numpy and take the per-site visit-count arrays
produced by a :class:`SignedWalker` instance (``walker.visits_agg``)
or by :func:`aggregate_visits`. Nothing in this module depends on the
walker class itself — an external analysis (notebook, post-processing
script) can call these helpers directly on any pair of length-``N``
non-negative integer arrays.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ...config.const import SRW_OVERLAP_KINDS


def aggregate_visits(visits_per_walker: np.ndarray) -> np.ndarray:
    """Sum per-walker visit counts across an ensemble.

    Parameters
    ----------
    visits_per_walker : (n_walkers, N) int array
        Per-walker, per-site visit counts.

    Returns
    -------
    visits_agg : (N,) int64
        Sum over walkers.
    """
    vw = np.asarray(visits_per_walker)
    if vw.ndim != 2:
        raise ValueError(
            f"expected (n_walkers, N) array, got shape {vw.shape}."
        )
    return vw.sum(axis=0, dtype=np.int64)


def _as_density(c: np.ndarray) -> np.ndarray:
    """Normalise a non-negative count array to a probability density.

    Returns a zero array if the input sums to zero.
    """
    c = np.asarray(c, dtype=np.float64)
    total = c.sum()
    if total <= 0:
        return np.zeros_like(c, dtype=np.float64)
    return c / total


def overlap(
    c1: np.ndarray,
    c2: np.ndarray,
    *,
    kind: str = 'l2_raw',
) -> float:
    r"""Scalar overlap between two visit-count fields.

    ``c1``, ``c2`` are length-``N`` non-negative arrays.

    Parameters
    ----------
    kind : str
        One of:

        * ``'l2_raw'`` (default) — :math:`\langle c_1, c_2\rangle =
          \sum_x c_1(x)\, c_2(x)`, the primary per-user request
          (raw per-site counts, no normalisation).
        * ``'l2_norm'`` — same but with each field normalised to a
          probability density (:math:`\rho_i = c_i / \sum c_i`).
        * ``'l2_rescaled'`` — :math:`N \cdot \langle \rho_1,\rho_2\rangle`,
          the IPR-like finite-size-scaling order parameter.
        * ``'fidelity'`` — :math:`\sum_x \sqrt{\rho_1(x)\, \rho_2(x)}`
          (Bhattacharyya coefficient, in ``[0, 1]``).
        * ``'tv'`` — :math:`\tfrac12\sum_x |\rho_1(x) - \rho_2(x)|`
          (total-variation distance, in ``[0, 1]``).
        * ``'jaccard'`` — binary Jaccard of the supports:
          :math:`|\{x: c_1>0\} \cap \{x: c_2>0\}| /
          |\{x: c_1>0\} \cup \{x: c_2>0\}|`.
    """
    if kind not in SRW_OVERLAP_KINDS:
        raise ValueError(
            f"kind must be one of {SRW_OVERLAP_KINDS}, got {kind!r}."
        )
    c1 = np.asarray(c1)
    c2 = np.asarray(c2)
    if c1.shape != c2.shape:
        raise ValueError(
            f"shape mismatch: {c1.shape} vs {c2.shape}."
        )
    if c1.ndim != 1:
        raise ValueError(
            f"expected 1-D arrays, got {c1.ndim}-D."
        )
    if (c1 < 0).any() or (c2 < 0).any():
        raise ValueError("visit counts must be non-negative.")
    if kind == 'l2_raw':
        return float(np.dot(c1.astype(np.float64), c2.astype(np.float64)))
    if kind == 'jaccard':
        s1 = c1 > 0
        s2 = c2 > 0
        inter = int((s1 & s2).sum())
        union = int((s1 | s2).sum())
        return float(inter / union) if union else 0.0
    rho1 = _as_density(c1)
    rho2 = _as_density(c2)
    if kind == 'l2_norm':
        return float(np.dot(rho1, rho2))
    if kind == 'l2_rescaled':
        return float(rho1.size * np.dot(rho1, rho2))
    if kind == 'fidelity':
        return float(np.sqrt(rho1 * rho2).sum())
    # kind == 'tv'
    return float(0.5 * np.abs(rho1 - rho2).sum())


def overlap_all(
    c1: np.ndarray,
    c2: np.ndarray,
) -> dict:
    """Return every supported overlap flavour as a single dict.

    Convenience helper for the CLI ``--mode pair`` path.
    """
    return {kind: overlap(c1, c2, kind=kind) for kind in SRW_OVERLAP_KINDS}


def pairwise_overlap_matrix(
    visits_per_walker: np.ndarray,
    *,
    kind: str = 'l2_raw',
) -> np.ndarray:
    """Gram-style matrix of overlaps across walkers in one ensemble.

    Parameters
    ----------
    visits_per_walker : (n_walkers, N) array
    kind : str — see :func:`overlap`.

    Returns
    -------
    (n_walkers, n_walkers) float64 array.
    """
    vw = np.asarray(visits_per_walker)
    n = vw.shape[0]
    out = np.empty((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i, n):
            o = overlap(vw[i], vw[j], kind=kind)
            out[i, j] = o
            out[j, i] = o
    return out


__all__ = [
    'aggregate_visits',
    'overlap',
    'overlap_all',
    'pairwise_overlap_matrix',
]
