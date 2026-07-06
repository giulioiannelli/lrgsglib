"""Model-level defaults for XYModel (legacy class + new scheme classes).

Single-sources every tunable of the XY dynamics layer (Hard-Ban #1: no magic
numbers at call sites). The axis vocabularies (rules, update modes, orders,
coupling norms) live in the shared spin-model homes (``statsys._thermal`` /
``statsys._csr``) and are re-exported here so XY code has one import home —
the same convention as ``PottsModel/defaults``.

Citations (``.agents/ref/00-REFERENCES.md`` §"MCMC / spin-dynamics kernels"):
the planar-rotator Hamiltonian ``H = −Σ J_ij cos(θ_i − θ_j)`` and the
magnetisation ``m = |Σ e^{iθ}|/N`` are the O(2) case of the collective-update
paper → M6 (Wolff 1989), whose embedded-Ising reflection also powers the XY
cluster moves; acceptance rules map to M1/M2 as in the Ising defaults.
"""

from __future__ import annotations

import numpy as np

from .._csr import COUPLING_NORMS
from .._thermal import (  # noqa: F401  (re-exported vocabulary)
    ASYNC_ORDERS,
    THERMAL_RULES,
    UPDATE_MODES,
)

# ---------------------------------------------------------------------------
# Solver registry
# ---------------------------------------------------------------------------
# Registry key under which XYModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``POTTS_SOLVER_NAME``). ONE key per
# model: the legacy ``XYModel`` and the new scheme classes (``XYBase``
# leaves) both dispatch through it; the solver polymorphs on the instance.
XY_SOLVER_NAME: str = "xy"

# ---------------------------------------------------------------------------
# Model axes (new scheme classes)
# ---------------------------------------------------------------------------
# Temperature (k_B = 1); the historical XYModel constructor default.
XY_T_DEFAULT: float = 1.0
# Maximum angular perturbation of the single-site proposal (radians): the
# proposal is uniform in [θ − delta, θ + delta] (mod 2π), symmetric, so every
# acceptance rule below satisfies detailed balance. Historical default.
XY_DELTA_DEFAULT: float = 0.5
# Acceptance rule: standard Metropolis (M1) preserves the historical kernel.
# 'heatbath' is not wired for a continuous state (the conditional resample
# from the local von-Mises distribution is not implemented); 'glauber'
# (Barker) stays valid for any symmetric proposal.
XY_RULE_DEFAULT: str = "metropolis"
# Update schedule: asynchronous single-site kinetics. 'gillespie' is not
# available for a continuous state (the shared BKL engine assumes a
# deterministic local proposal; a continuous proposal is a fresh random draw).
XY_UPD_MODE_DEFAULT: str = "async"
# Async site-visit order: uniform random WITH replacement (true asynchronous
# kinetics). The legacy Python loop used a fresh permutation per sweep —
# kept as an explicit option, no longer the default (plan §3b).
XY_ORDER_DEFAULT: str = "random"
# Elementary moves: single-angle perturbation; Wolff (M6) / Swendsen–Wang
# (M5) cluster moves via M6's embedded-Ising construction — reflect the
# cluster's spins about a random direction. A reflection is an involution
# (like the Ising flip and unlike the Potts relabel), so the embedded bond
# satisfaction is flip-covariant and cluster moves stay valid on SIGNED
# couplings. Kawasaki pair exchange (M8) conserves the multiset of angles.
# No 'spectral' move: the eigenvector-subspace walk is an Ising (±1)
# construction.
XY_MOVES: tuple[str, ...] = ("single", "wolff", "sw", "exchange")
XY_MOVE_DEFAULT: str = "single"

# ---------------------------------------------------------------------------
# Hamiltonian conventions (plan §3b, shared coupling builder statsys._csr)
# ---------------------------------------------------------------------------
XY_COUPLING_NORMS: tuple[str, ...] = COUPLING_NORMS
XY_COUPLING_NORM_DEFAULT: str = "raw"

# ---------------------------------------------------------------------------
# Observables (new scheme classes; ObservableSet names + file bases)
# ---------------------------------------------------------------------------
XY_OBS_ENERGY: str = "energy"  # intensive per-site energy series
XY_OBS_MAGN: str = "magn"  # magnetisation series |Σ e^{iθ}|/N
XY_OBS_SNAPSHOTS: str = "snapshots"  # (n_rec, N) float64 angle trajectory

# Default recorded set — the legacy sampler recorded energy + magnetisation
# when save_observables was on; snapshots are opt-in (I/O-heavy).
XY_OBS_DEFAULT: tuple[str, ...] = (XY_OBS_ENERGY, XY_OBS_MAGN)

# On-disk file bases (``<base>_p=<pflip>[_suffix].bin`` under dynpath),
# matching the Ising/Potts token conventions. dynpath is
# ``<path_data>/XY_DYN_SUBDIR`` (the legacy XYModel location).
XY_ENERGY_FBASE: str = "ene"
XY_MAGN_FBASE: str = "m"
XY_SNAPSHOTS_FBASE: str = "sout"
XY_DYN_SUBDIR: str = "xy"

# Snapshot cadence + size guard (same shape as Ising/Potts): keep every k-th
# sweep when streaming, refuse before writing past the byte cap (float64
# angles = 8 bytes/site/record).
XY_SNAPSHOT_EVERY_DEFAULT: int = 1
XY_SNAPSHOT_MAX_BYTES: int = 512 * 1024 * 1024  # 512 MiB, as Ising/Potts

# Full circle (radians): the angle state lives in [0, XY_TWO_PI).
XY_TWO_PI: float = 2.0 * np.pi
