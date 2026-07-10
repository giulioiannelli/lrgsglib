"""Default constants for ContactProcess dynamics + its disk-backed observables.

Mirrors ``VoterModel/defaults.py`` so the contact process exposes the same
ObservableSet-driven persistence surface (Phase 2/4 of the modernization spine):
a ``density`` scalar series (active-site fraction over MCMC time) and a
``snapshots`` spin trajectory, either kept in RAM (``savedisk=False``) or
streamed/dumped under ``dynpath`` and exposed as lazy disk-backed attributes.
"""

from __future__ import annotations

__all__ = [
    "CP_SOLVER_NAME",
    "CP_OBS_DENSITY",
    "CP_OBS_SNAPSHOTS",
    "CP_DENSITY_FBASE",
    "CP_SNAPSHOTS_FBASE",
    "CP_SNAPSHOT_EVERY",
    "CP_SNAPSHOT_MAX_BYTES",
    "CP_SIR_BETA_DEFAULT",
    "CP_SIR_MU_DEFAULT",
    "CP_EI_ACTIVATION_DEFAULT",
    "CP_EI_NLS_DEFAULT",
    "CP_EI_UPDATE_MODE_DEFAULT",
]

# ----------------------------------------------------------------------------
# Rate defaults for the SIR (infection/recovery) contact process. ``beta`` is
# the infection rate accumulated per active positive-edge neighbour; ``mu`` is
# the spontaneous recovery rate. beta=mu=1 reproduces the historical mu-only
# behaviour where the infection rate equalled the edge weight.
CP_SIR_BETA_DEFAULT: float = 1.0
CP_SIR_MU_DEFAULT: float = 1.0

# ----------------------------------------------------------------------------
# EI (excitation-inhibition) defaults: the C kernels' non-linearity, the
# number of log-spaced density/snapshot samples, and the update strategy.
# These are also the elision defaults of the ``act=`` / ``nls=`` / ``upd=``
# run-dirname tokens (Phase-C naming).
CP_EI_ACTIVATION_DEFAULT: str = "tanh"
CP_EI_NLS_DEFAULT: int = 1000
CP_EI_UPDATE_MODE_DEFAULT: str = "naive"

# Solver-registry name for ContactProcess (key in ``statsys._solver_engine``).
CP_SOLVER_NAME: str = "contact_process"

# ObservableSet keys (the names used to index ``self.observables``).
CP_OBS_DENSITY: str = "density"
CP_OBS_SNAPSHOTS: str = "snapshots"

# On-disk filename bases (joined with the run's out_suffix by ``get_p_fname``).
CP_DENSITY_FBASE: str = "rho"  # active-site density series (float64 .bin)
CP_SNAPSHOTS_FBASE: str = (
    "sout"  # spin-configuration trajectory (int8 .bin, (n_rec, N))
)

# Default snapshot subsampling: keep every k-th recorded sweep (1 == every sweep).
CP_SNAPSHOT_EVERY: int = 1

# Hard guard on the (subsampled) snapshot trajectory size; the run raises BEFORE
# writing if exceeded (raise the cap or snapshot_every). Mirrors VoterModel.
CP_SNAPSHOT_MAX_BYTES: int = 512 * 1024 * 1024  # 512 MiB
