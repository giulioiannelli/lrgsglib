"""CLI defaults for the L2D_VoterModel program and its serializer.

The voter *dynamics* tunables (rule, q, eps, alpha, upd_mode, cluster_mode,
absorbing_every, sout_every, ...) are single-sourced in
``lrgsglib.statsys.VoterModel.defaults`` and imported by ``progargs.VoterModel``
(project hard-ban #1: no duplicated tunables). The generic dispatch / array /
quench / save-results defaults are reused from ``defs.IsingDynamics``. Only the
program/serializer CLI knobs with no canonical home elsewhere live here.
"""

from __future__ import annotations

# --- program backend / init / time (CLI-layer defaults) ---
DEFAULT_VOTER_RUNLANG: str = "py"
DEFAULT_VOTER_INIT_COND: str = "uniform"
DEFAULT_VOTER_STATE_TYPE: str = "bipolar"
DEFAULT_VOTER_STEPS: int | None = None
DEFAULT_VOTER_SIMREF: float | None = None
DEFAULT_VOTER_FREQ: int = 10
DEFAULT_VOTER_NUM_LOG_SAMPLES: int = 100
DEFAULT_VOTER_SEED: int | None = None
DEFAULT_VOTER_SOUT_NLOG: int | None = None
DEFAULT_VOTER_SOUT_FORCE_FULL: bool = False
DEFAULT_VOTER_RANDSTR: bool = True

# --- observable save gates (program defaults) ---
# Magnetization is the headline observable -> on by default; the trajectory
# (sout) and cluster-size distribution (cldist) are heavier -> opt-in.
DEFAULT_VOTER_SAVE_MAGN: bool = True
DEFAULT_VOTER_SAVE_SOUT: bool = False
DEFAULT_VOTER_SAVE_CLDIST: bool = False
DEFAULT_VOTER_SAVE_ALL: bool = False
