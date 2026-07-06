"""Model-level defaults for HeisenbergModel (legacy class + new scheme classes).

Single-sources every tunable of the Heisenberg dynamics layer (Hard-Ban #1:
no magic numbers at call sites). The axis vocabularies (rules, update modes,
orders, coupling norms) live in the shared spin-model homes
(``statsys._thermal`` / ``statsys._csr``) and are re-exported here so
Heisenberg code has one import home — the same convention as
``XYModel/defaults``.

Citations (``.agents/ref/00-REFERENCES.md`` §"MCMC / spin-dynamics kernels"):
the classical-Heisenberg Hamiltonian ``H = −Σ J_ij (n_i · n_j)`` and the
magnetisation ``m = |Σ n_i|/N`` are the O(3) case of the collective-update
paper → M6 (Wolff 1989), whose embedded-Ising reflection also powers the
cluster moves; acceptance rules map to M1/M2 as in the Ising defaults.
"""

from __future__ import annotations

from .._csr import COUPLING_NORMS
from .._thermal import (  # noqa: F401  (re-exported vocabulary)
    ASYNC_ORDERS,
    THERMAL_RULES,
    UPDATE_MODES,
)

# ---------------------------------------------------------------------------
# Solver registry
# ---------------------------------------------------------------------------
# Registry key under which HeisenbergModel's solver backends are registered
# in ``statsys._solver_engine`` (mirrors ``XY_SOLVER_NAME``). ONE key per
# model: the legacy ``HeisenbergModel`` and the new scheme classes
# (``HeisenbergBase`` leaves) both dispatch through it; the solver polymorphs
# on the instance.
HEISENBERG_SOLVER_NAME: str = "heisenberg"

# ---------------------------------------------------------------------------
# Model axes (new scheme classes)
# ---------------------------------------------------------------------------
# Temperature (k_B = 1); the historical HeisenbergModel constructor default.
HEISENBERG_T_DEFAULT: float = 1.0
# Perturbation amplitude of the single-site proposal: the proposal is
# ``normalize(n + delta·u)`` with ``u`` uniform in [−1, 1]³ — the historical
# kernel shared by the legacy Python loop and the native/C kernels.
HEISENBERG_DELTA_DEFAULT: float = 0.3
# Acceptance rule: standard Metropolis (M1) preserves the historical kernel.
# 'heatbath' is not wired for a continuous state; 'glauber' (Barker) stays
# valid for the same proposal the metropolis rule uses.
HEISENBERG_RULE_DEFAULT: str = "metropolis"
# Update schedule: asynchronous single-site kinetics. 'gillespie' is not
# available for a continuous state (the shared BKL engine assumes a
# deterministic local proposal).
HEISENBERG_UPD_MODE_DEFAULT: str = "async"
# Async site-visit order: uniform random WITH replacement (true asynchronous
# kinetics). The legacy Python loop used a fresh permutation per sweep —
# kept as an explicit option, no longer the default (plan §3b).
HEISENBERG_ORDER_DEFAULT: str = "random"
# Elementary moves: single-spin perturbation; Wolff (M6) / Swendsen–Wang
# (M5) cluster moves via M6's embedded-Ising construction — reflect the
# cluster's spins about a random unit vector. A reflection is an involution
# (like the Ising flip), so the embedded bond satisfaction is flip-covariant
# and cluster moves stay valid on SIGNED couplings. Kawasaki pair exchange
# (M8) conserves the multiset of spin vectors. No 'spectral' move: the
# eigenvector-subspace walk is an Ising (±1) construction.
HEISENBERG_MOVES: tuple[str, ...] = ("single", "wolff", "sw", "exchange")
HEISENBERG_MOVE_DEFAULT: str = "single"

# ---------------------------------------------------------------------------
# Hamiltonian conventions (plan §3b, shared coupling builder statsys._csr)
# ---------------------------------------------------------------------------
HEISENBERG_COUPLING_NORMS: tuple[str, ...] = COUPLING_NORMS
HEISENBERG_COUPLING_NORM_DEFAULT: str = "raw"

# ---------------------------------------------------------------------------
# Observables (new scheme classes; ObservableSet names + file bases)
# ---------------------------------------------------------------------------
HEISENBERG_OBS_ENERGY: str = "energy"  # intensive per-site energy series
HEISENBERG_OBS_MAGN: str = "magn"  # magnetisation series |Σ n_i|/N
HEISENBERG_OBS_SNAPSHOTS: str = "snapshots"  # (n_rec, 3N) float64 trajectory

# Default recorded set — the legacy sampler recorded energy + magnetisation
# when save_observables was on; snapshots are opt-in (I/O-heavy).
HEISENBERG_OBS_DEFAULT: tuple[str, ...] = (
    HEISENBERG_OBS_ENERGY,
    HEISENBERG_OBS_MAGN,
)

# On-disk file bases (``<base>_p=<pflip>[_suffix].bin`` under dynpath),
# matching the Ising/Potts/XY token conventions. dynpath is
# ``<path_data>/HEISENBERG_DYN_SUBDIR`` (the legacy HeisenbergModel location).
HEISENBERG_ENERGY_FBASE: str = "ene"
HEISENBERG_MAGN_FBASE: str = "m"
HEISENBERG_SNAPSHOTS_FBASE: str = "sout"
HEISENBERG_DYN_SUBDIR: str = "heisenberg"

# Snapshot cadence + size guard (same shape as Ising/Potts/XY): keep every
# k-th sweep when streaming, refuse before writing past the byte cap
# (float64 3-vectors = 24 bytes/site/record; a snapshot row is the flattened
# (N, 3) state, 3N doubles).
HEISENBERG_SNAPSHOT_EVERY_DEFAULT: int = 1
HEISENBERG_SNAPSHOT_MAX_BYTES: int = 512 * 1024 * 1024  # 512 MiB
