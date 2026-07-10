"""Vectorized sublattice-parallel (sync) sweeps — the np/cu backends.

Same M4 dynamics as ``ThermalEngine._compile_sync`` (greedy color
classes; per class: propose + ΔE from the frozen state, then independent
accept/commit — no intra-class bonds, so the block preserves the Gibbs
measure), but each color class is processed as ONE array operation.
True async (random-sequential) has a serial dependency per flip and is
NOT offered here; ``upd_mode='sync'`` is the vectorizable schedule, and
on bipartite lattices the color classes ARE the checkerboard.

Written against a generic array module ``xp`` so a single implementation
serves both backends: ``xp=numpy`` (np) and ``xp=cupy`` (cu).

A consuming substrate provides three vectorized hooks:

- ``_vec_propose(cls, xp, rng)`` -> proposals for the class sites;
- ``_vec_delta_E(cls, proposals, xp)`` -> ΔE array (frozen state);
- ``_vec_commit(cls, proposals, accept)`` -> apply accepted proposals,
  return the summed committed ΔE is handled here from the ΔE array.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

__all__ = ["greedy_color_classes", "VectorSyncEngine"]


def greedy_color_classes(
    nbr_ptr: np.ndarray, nbr_idx: np.ndarray, n_sites: int
) -> list[np.ndarray]:
    """Greedy-color the interaction graph: no two sites in a class are
    adjacent (same algorithm/order as ``ThermalEngine``)."""
    colors = np.full(n_sites, -1, dtype=np.int64)
    for i in range(n_sites):
        nbrs = nbr_idx[nbr_ptr[i] : nbr_ptr[i + 1]]
        used = {int(c) for c in colors[nbrs] if c >= 0}
        c = 0
        while c in used:
            c += 1
        colors[i] = c
    return [np.flatnonzero(colors == c) for c in range(int(colors.max()) + 1)]


class VectorSyncEngine:
    """Array-parallel synchronous kinetics over a vector-capable substrate.

    Parameters mirror ``ThermalEngine`` where they apply: ``rule`` is
    'metropolis' or 'glauber' (heat-bath resampling is not a
    propose/accept scheme and stays python-only), ``tie_flip_p`` is the
    discrete-substrate tie policy (ΔE == 0 under metropolis).
    """

    def __init__(
        self,
        substrate: Any,
        *,
        rule: str,
        T: float,
        tie_flip_p: float | None = None,
        xp: Any = np,
        rng: Any = None,
    ) -> None:
        if rule not in ("metropolis", "glauber"):
            raise NotImplementedError(
                f"VectorSyncEngine implements 'metropolis' and 'glauber'; "
                f"rule={rule!r} is python-only."
            )
        if T < 0.0:
            raise ValueError(f"Temperature must be >= 0; got {T}.")
        self.substrate = substrate
        self.rule = rule
        self.T = float(T)
        self.tie_flip_p = 1.0 if tie_flip_p is None else float(tie_flip_p)
        self.xp = xp
        self.rng = rng if rng is not None else xp.random

    def _accept_mask(self, dE: Any) -> Any:
        xp, T = self.xp, self.T
        u = self.rng.random(dE.shape[0])
        if self.rule == "glauber":
            if T > 0.0:
                p = 1.0 / (1.0 + xp.exp(dE / T))
            else:
                p = xp.where(dE < 0.0, 1.0, xp.where(dE > 0.0, 0.0, 0.5))
            return u < p
        # metropolis
        if T > 0.0:
            p = xp.exp(-xp.maximum(dE, 0.0) / T)
            accept = u < p
        else:
            accept = dE < 0.0
        ties = dE == 0.0
        if bool(xp.any(ties)):
            accept = xp.where(
                ties,
                self.rng.random(dE.shape[0]) < self.tie_flip_p,
                accept,
            )
        return accept

    def compile_step(self) -> Callable[[], float]:
        substrate = self.substrate
        xp = self.xp
        classes = [
            xp.asarray(cls)
            for cls in greedy_color_classes(
                substrate._nbr_ptr, substrate._nbr_idx, substrate.N
            )
        ]

        def sweep() -> float:
            dE_sweep = 0.0
            for cls in classes:
                proposals = substrate._vec_propose(cls, xp, self.rng)
                dE = substrate._vec_delta_E(cls, proposals, xp)
                accept = self._accept_mask(dE)
                substrate._vec_commit(cls, proposals, accept)
                dE_sweep += float(dE[accept].sum())
            return dE_sweep

        return sweep
