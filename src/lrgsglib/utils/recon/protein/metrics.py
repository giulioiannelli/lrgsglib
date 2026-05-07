"""Lightweight protein-reconstruction metrics: Kabsch RMSD and Q3."""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

__all__ = ["kabsch_align", "kabsch_rmsd", "q3_score"]


def kabsch_align(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return ``P`` rotated and translated onto ``Q`` (least-squares CA superposition)."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return Pc @ R.T + Q.mean(axis=0)


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Optimal-superposition RMSD between two ``(n, 3)`` coordinate sets."""
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    aligned = kabsch_align(P, Q)
    diff = aligned - Q
    return float(np.sqrt((diff * diff).sum() / max(len(P), 1)))


def q3_score(ss_pred: Sequence[str] | Iterable[str],
             ss_true: Sequence[str] | Iterable[str]) -> float:
    """Fraction of residues where predicted SS class matches ground truth.

    Inputs are sequences of single-letter SS codes (``H``, ``E``, ``C`` per
    the existing ``assign_secondary_structure_from_coords`` convention).
    """
    pred = list(ss_pred)
    true = list(ss_true)
    n = min(len(pred), len(true))
    if n == 0:
        return 0.0
    return sum(p == t for p, t in zip(pred[:n], true[:n])) / n
