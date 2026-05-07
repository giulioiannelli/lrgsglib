"""Coordinate-padding utility used by protein TMD encoders.

The previous implementation hosted three near-duplicate basis builders; those
have been collapsed into :mod:`substrate` (see ``build_substrate_basis``).
Only the lightweight ``pad_protein_coordinates`` helper survives here, since
it is consumed by the encode lab and by tests.
"""

from __future__ import annotations

import numpy as np

__all__ = ["pad_protein_coordinates"]


def pad_protein_coordinates(coords: np.ndarray, target_residues: int) -> np.ndarray:
    """Pad or truncate an ``(n, 3)`` coordinate array to exactly ``target_residues`` rows.

    Padding rows are filled with the centroid of ``coords`` so that the
    centred-coordinate signal projected on the constant Laplacian eigenmode
    remains numerically zero.
    """
    n_residues = coords.shape[0]
    if n_residues == target_residues:
        return coords
    if n_residues < target_residues:
        padded_coords = np.zeros((target_residues, 3), dtype=coords.dtype)
        padded_coords[:n_residues] = coords
        if n_residues > 0:
            centroid = coords.mean(axis=0)
            padded_coords[n_residues:] = centroid
        return padded_coords
    return coords[:target_residues]
