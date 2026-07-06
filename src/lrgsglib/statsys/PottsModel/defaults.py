"""Model-level defaults for PottsModel (legacy class + new scheme classes).

Single-sources every tunable of the Potts dynamics layer (Hard-Ban #1: no
magic numbers at call sites). The axis vocabularies (rules, update modes,
orders, coupling norms) live in the shared spin-model homes
(``statsys._thermal`` / ``statsys._csr``) and are re-exported here so Potts
code has one import home — the same convention as ``IsingDynamics/defaults``.

Citations (``.agents/ref/00-REFERENCES.md`` §"MCMC / spin-dynamics kernels"):
the q-state Potts Hamiltonian ``H = −Σ J_ij δ(σ_i, σ_j)`` and the order
parameter ``(q·max_frac − 1)/(q − 1)`` follow the Potts-model review → M10
(Wu 1982); acceptance rules and cluster moves map to M1/M2/M5/M6/M8 as in the
Ising defaults.
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
# Registry key under which PottsModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``VOTER_SOLVER_NAME``). ONE key per
# model: the legacy ``PottsModel`` and the new scheme classes (``PottsBase``
# leaves) both dispatch through it; the solver polymorphs on the instance.
POTTS_SOLVER_NAME: str = "potts"

# ---------------------------------------------------------------------------
# Model axes (new scheme classes)
# ---------------------------------------------------------------------------
# Number of Potts states; q = 2 with the bipolar mapping {0,1} -> {-1,+1} is
# the Ising model at halved couplings (δ(σi,σj) = (1 + s_i s_j)/2, M10).
POTTS_Q_DEFAULT: int = 3
# Temperature (k_B = 1); the historical PottsModel constructor default.
POTTS_T_DEFAULT: float = 1.0
# Acceptance rule: standard Metropolis (M1) preserves the historical kernel.
# 'heatbath' is q = 2 only (the shared kernel is the two-state conditional
# resample; the multi-state conditional is not wired).
POTTS_RULE_DEFAULT: str = "metropolis"
# Update schedule: asynchronous single-site kinetics. 'gillespie' is q = 2
# only (the shared BKL engine assumes a deterministic local proposal; q > 2
# needs per-(site,label) rate channels).
POTTS_UPD_MODE_DEFAULT: str = "async"
# Async site-visit order: uniform random WITH replacement (true asynchronous
# kinetics). The legacy Python loop used a fresh permutation per sweep —
# kept as an explicit option, no longer the default (plan §3b).
POTTS_ORDER_DEFAULT: str = "random"
# Elementary moves: single-label change; Wolff (M6) / Swendsen–Wang (M5)
# cluster relabelings — non-negative couplings only (a joint relabel does not
# preserve an unlike pair, so antiferromagnetic bonds break the FK
# bookkeeping, unlike the Ising joint flip); Kawasaki pair exchange (M8),
# conserving the per-state densities. No 'spectral' move: the eigenvector-
# subspace walk is an Ising (±1) construction.
POTTS_MOVES: tuple[str, ...] = ("single", "wolff", "sw", "exchange")
POTTS_MOVE_DEFAULT: str = "single"
# Probability of accepting a tie (ΔE == 0) under rule='metropolis' (same
# project kinetic extension as Ising; presets in TIE_FLIP_PRESETS). 1.0 is
# the standard M1 behaviour and matches the native kernel, which accepts
# every ΔE <= 0 proposal at ALL temperatures (including T = 0).
POTTS_TIE_FLIP_DEFAULT: float = TIE_FLIP_PRESETS["standard"]

# ---------------------------------------------------------------------------
# Hamiltonian conventions (plan §3b, shared coupling builder statsys._csr)
# ---------------------------------------------------------------------------
POTTS_COUPLING_NORMS: tuple[str, ...] = COUPLING_NORMS
POTTS_COUPLING_NORM_DEFAULT: str = "raw"

# ---------------------------------------------------------------------------
# Observables (new scheme classes; ObservableSet names + file bases)
# ---------------------------------------------------------------------------
POTTS_OBS_ENERGY: str = "energy"  # intensive per-site energy series
POTTS_OBS_MAGN: str = "magn"  # order-parameter series (q·max_frac − 1)/(q − 1)
POTTS_OBS_SNAPSHOTS: str = "snapshots"  # (n_rec, N) int32 label trajectory

# Default recorded set — the legacy sampler recorded energy + order parameter
# when save_observables was on; snapshots are opt-in (I/O-heavy).
POTTS_OBS_DEFAULT: tuple[str, ...] = (POTTS_OBS_ENERGY, POTTS_OBS_MAGN)

# On-disk file bases (``<base>_p=<pflip>[_suffix].bin`` under dynpath),
# matching the Ising/Voter/CP token conventions. dynpath is
# ``<path_data>/POTTS_DYN_SUBDIR`` (the legacy PottsModel location).
POTTS_ENERGY_FBASE: str = "ene"
POTTS_MAGN_FBASE: str = "m"
POTTS_SNAPSHOTS_FBASE: str = "sout"
POTTS_DYN_SUBDIR: str = "potts"

# Snapshot cadence + size guard (same shape as Ising/CP): keep every k-th
# sweep when streaming, refuse before writing past the byte cap (int32
# labels = 4 bytes/site/record).
POTTS_SNAPSHOT_EVERY_DEFAULT: int = 1
POTTS_SNAPSHOT_MAX_BYTES: int = 512 * 1024 * 1024  # 512 MiB, as Ising/CP
