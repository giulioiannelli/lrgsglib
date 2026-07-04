"""Spectral proposal / field construction shared by dynamics models.

PHASE-0 SKELETON — shapes only. The bodies are extracted from today's Ising
monolith (the ``topo_*`` runlang paths) during the Ising vertical slice
(Phase 1), dtype-agnostic so the same helpers serve ``BinDynSys``- and
``VecDynSys``-based models.

Two capabilities live here (plan §4):

- **eigenvector-subspace proposal** — propose moves inside the span of the
  low-lying Laplacian eigenvectors (today's "topological" Metropolis / CEM
  moves, ``move='spectral'`` after the refactor);
- **RBIM field construction** — build the external field from a (binarized)
  eigenvector of the signed Laplacian (today's TFCA field,
  ``field='spectral'`` on the annealing scheme).

Both consume the graph through ``sg.get_eigV_bin_check`` /
``sg.compute_laplacian_spectrum_weigV`` (see
``graphs.protocols.DynamicsGraphProtocol`` / ``SpectralGraphProtocol``); no
model internals leak in here.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import DTypeLike, NDArray

if TYPE_CHECKING:
    from ..graphs.protocols import DynamicsGraphProtocol

__all__ = ["eigvec_subspace_proposal", "build_spectral_field"]


def eigvec_subspace_proposal(
    sg: "DynamicsGraphProtocol",
    n_modes: int,
    dtype: DTypeLike = np.float64,
) -> NDArray:
    """Build the proposal basis spanned by the *n_modes* lowest Laplacian
    eigenvectors (Phase-1 body; extracted from the Ising ``topo_*`` paths)."""
    raise NotImplementedError(
        "Phase-0 skeleton: extracted from the Ising monolith in Phase 1."
    )


def build_spectral_field(
    sg: "DynamicsGraphProtocol",
    which: int = 0,
    amplitude: float = 1.0,
    dtype: DTypeLike = np.float64,
) -> NDArray:
    """Build an external field aligned with eigenvector *which* of the signed
    Laplacian (today's TFCA field; Phase-1 body)."""
    raise NotImplementedError(
        "Phase-0 skeleton: extracted from the Ising monolith in Phase 1."
    )
