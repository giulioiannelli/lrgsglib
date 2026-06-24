"""Default constants for ContactProcess dynamics + its disk-backed observables.

Mirrors ``VoterModel/defaults.py`` so the contact process exposes the same
ObservableSet-driven persistence surface (Phase 2/4 of the modernization spine):
a ``density`` scalar series (active-site fraction over MCMC time) and a
``snapshots`` spin trajectory, either kept in RAM (``savedisk=False``) or
streamed/dumped under ``dynpath`` and exposed as lazy disk-backed attributes.
"""

from __future__ import annotations

__all__ = [
    "CP_OBS_DENSITY",
    "CP_OBS_SNAPSHOTS",
    "CP_DENSITY_FBASE",
    "CP_SNAPSHOTS_FBASE",
    "CP_SNAPSHOT_EVERY",
    "CP_SNAPSHOT_MAX_BYTES",
]

# ObservableSet keys (the names used to index ``self.observables``).
CP_OBS_DENSITY: str = "density"
CP_OBS_SNAPSHOTS: str = "snapshots"

# On-disk filename bases (joined with the run's out_suffix by ``get_p_fname``).
CP_DENSITY_FBASE: str = "rho"      # active-site density series (float64 .bin)
CP_SNAPSHOTS_FBASE: str = "sout"   # spin-configuration trajectory (int8 .bin, (n_rec, N))

# Default snapshot subsampling: keep every k-th recorded sweep (1 == every sweep).
CP_SNAPSHOT_EVERY: int = 1

# Hard guard on the (subsampled) snapshot trajectory size; the run raises BEFORE
# writing if exceeded (raise the cap or snapshot_every). Mirrors VoterModel.
CP_SNAPSHOT_MAX_BYTES: int = 512 * 1024 * 1024     # 512 MiB
