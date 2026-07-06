"""Model-level defaults for MultiSpeciesModel (legacy + new scheme classes).

Single-sources every tunable of the multi-species dynamics layer (Hard-Ban
#1: no magic numbers at call sites). The axis vocabularies (rules, update
modes, orders, coupling norms) live in the shared spin-model homes
(``statsys._thermal`` / ``statsys._csr``) and are re-exported here so
multi-species code has one import home — the same convention as
``PottsModel/defaults``.

Citations (``.agents/ref/00-REFERENCES.md`` §"MCMC / spin-dynamics kernels"):
the model is k coupled q-state Potts components per site,
``H = −Σ_{(i,j)} Σ_{s1,s2} J_ij M[s1,s2] δ(σ_i^{s1}, σ_j^{s2})`` — the
single-species case (k = 1, M = 1) is exactly the Potts model → M10 (Wu
1982); acceptance rules map to M1/M2 as in the Potts defaults.
"""

from __future__ import annotations

from .._csr import COUPLING_NORMS
from .._thermal import (  # noqa: F401  (re-exported vocabulary)
    ASYNC_ORDERS,
    THERMAL_RULES,
    TIE_FLIP_PRESETS,
    UPDATE_MODES,
)

# ---------------------------------------------------------------------------
# Solver registry
# ---------------------------------------------------------------------------
# Registry key under which MultiSpeciesModel's solver backends are registered
# in ``statsys._solver_engine`` (mirrors ``POTTS_SOLVER_NAME``). ONE key per
# model: the legacy ``MultiSpeciesModel`` and the new scheme classes
# (``MultiSpeciesBase`` leaves) both dispatch through it; the solver
# polymorphs on the instance.
MULTISPEC_SOLVER_NAME: str = "multispec"

# ---------------------------------------------------------------------------
# Model axes (new scheme classes)
# ---------------------------------------------------------------------------
# Number of coupled species components per site and states per component;
# the historical MultiSpeciesModel constructor defaults. k = 1 with M = 1 is
# exactly the q-state Potts model (M10).
MULTISPEC_SPECIES_DEFAULT: int = 2
MULTISPEC_Q_DEFAULT: int = 3
# Temperature (k_B = 1); the historical constructor default.
MULTISPEC_T_DEFAULT: float = 1.0
# Acceptance rule: standard Metropolis (M1) preserves the historical kernel.
# 'heatbath' is not wired (the conditional resample over the coupled
# (component, label) channels is not implemented); 'glauber' (Barker) stays
# valid for the symmetric (component, label) proposal.
MULTISPEC_RULE_DEFAULT: str = "metropolis"
# Update schedule: asynchronous single-site kinetics. 'gillespie' is not
# available (the shared BKL engine assumes a deterministic local proposal;
# the multi-species proposal draws a random component AND label).
MULTISPEC_UPD_MODE_DEFAULT: str = "async"
# Async site-visit order: uniform random WITH replacement (true asynchronous
# kinetics). The legacy Python loop used a fresh permutation per sweep —
# kept as an explicit option, no longer the default (plan §3b).
MULTISPEC_ORDER_DEFAULT: str = "random"
# Elementary moves: single-component relabel; Kawasaki full-site swap (M8),
# conserving every per-species density. No cluster moves: the FK
# construction for the COUPLED multi-species Hamiltonian (inter-species
# delta terms) is not wired. No 'spectral' move (Ising ±1 construction).
MULTISPEC_MOVES: tuple[str, ...] = ("single", "exchange")
MULTISPEC_MOVE_DEFAULT: str = "single"
# Probability of accepting a tie (ΔE == 0) under rule='metropolis' (same
# project kinetic extension as Ising/Potts; presets in TIE_FLIP_PRESETS).
# 1.0 is the standard M1 behaviour and matches the native kernel, which
# accepts every ΔE <= 0 proposal at ALL temperatures (including T = 0).
MULTISPEC_TIE_FLIP_DEFAULT: float = TIE_FLIP_PRESETS["standard"]

# ---------------------------------------------------------------------------
# Hamiltonian conventions (plan §3b, shared coupling builder statsys._csr)
# ---------------------------------------------------------------------------
MULTISPEC_COUPLING_NORMS: tuple[str, ...] = COUPLING_NORMS
MULTISPEC_COUPLING_NORM_DEFAULT: str = "raw"

# ---------------------------------------------------------------------------
# Observables (new scheme classes; ObservableSet names + file bases)
# ---------------------------------------------------------------------------
MULTISPEC_OBS_ENERGY: str = "energy"  # intensive per-site energy series
MULTISPEC_OBS_MAGN: str = "magn"  # species-averaged Potts order parameter
MULTISPEC_OBS_SNAPSHOTS: str = "snapshots"  # (n_rec, kN) int32 trajectory

# Default recorded set — the legacy sampler recorded energy only; the
# species-averaged order parameter joins the cross-model convention.
# Snapshots are opt-in (I/O-heavy).
MULTISPEC_OBS_DEFAULT: tuple[str, ...] = (
    MULTISPEC_OBS_ENERGY,
    MULTISPEC_OBS_MAGN,
)

# On-disk file bases (``<base>_p=<pflip>[_suffix].bin`` under dynpath),
# matching the Ising/Potts token conventions. dynpath is
# ``<path_data>/MULTISPEC_DYN_SUBDIR`` (the legacy location).
MULTISPEC_ENERGY_FBASE: str = "ene"
MULTISPEC_MAGN_FBASE: str = "m"
MULTISPEC_SNAPSHOTS_FBASE: str = "sout"
MULTISPEC_DYN_SUBDIR: str = "multi_species"

# Snapshot cadence + size guard (same shape as the other spin models): keep
# every k-th sweep when streaming, refuse before writing past the byte cap
# (int32 labels = 4 bytes/component/record; a snapshot row is the flattened
# (N, species) state).
MULTISPEC_SNAPSHOT_EVERY_DEFAULT: int = 1
MULTISPEC_SNAPSHOT_MAX_BYTES: int = 512 * 1024 * 1024  # 512 MiB
