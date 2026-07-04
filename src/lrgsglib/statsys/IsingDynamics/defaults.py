"""Model-level defaults for IsingDynamics (legacy class + new scheme classes).

Single-sources every tunable of the Ising dynamics layer (Hard-Ban #1: no
magic numbers at call sites). The axis vocabularies (rules, update modes,
orders) live in ``statsys._thermal`` — the shared spin-model engine — and are
re-exported here so Ising code has one import home.

Citations: acceptance rules and tie presets map to
``.agents/ref/00-REFERENCES.md`` §"MCMC / spin-dynamics kernels" (M1
Metropolis 1953, M2 Glauber 1963, M4 Newman & Barkema 1999).
"""

from __future__ import annotations

from .._thermal import (  # noqa: F401  (re-exported vocabulary)
    ASYNC_ORDERS,
    THERMAL_RULES,
    TIE_FLIP_PRESETS,
    UPDATE_MODES,
)

# ---------------------------------------------------------------------------
# Solver registry
# ---------------------------------------------------------------------------
# Registry key under which IsingDynamics' solver backends (py/pb/cu/C) are
# registered in ``statsys._solver_engine``. ONE key per model: the scheme
# classes (Metropolis/SA/PT/CEM) are the polymorphism the solver's
# ``execute(model)`` rides on.
ISING_SOLVER_NAME: str = "ising"

# ---------------------------------------------------------------------------
# Layer-2 axis defaults (plan §3b / §9)
# ---------------------------------------------------------------------------
# Acceptance rule: standard Metropolis (M1) preserves the historical kernel.
ISING_RULE_DEFAULT: str = "metropolis"
# Update schedule: asynchronous single-site kinetics.
ISING_UPD_MODE_DEFAULT: str = "async"
# Async site-visit order: uniform random WITH replacement (true asynchronous
# kinetics). The legacy Python loop was 'typewriter' (fixed 0..N-1) — kept as
# an explicit option, no longer the default (plan §3b correctness fix).
ISING_ORDER_DEFAULT: str = "random"
# Elementary move vocabulary; Phase 1 wires 'single' only, the rest are
# validated at setup and land with their phases (wolff/sw = cluster,
# spectral = eigvec-subspace via statsys._spectral, exchange = Kawasaki).
ISING_MOVES: tuple[str, ...] = ("single", "wolff", "sw", "spectral", "exchange")
ISING_MOVE_DEFAULT: str = "single"

# Probability of accepting a tie (ΔE == 0) spin flip under rule='metropolis'.
# 1.0 is standard Metropolis (exp(-0/T) = 1, M1); 0.5 is the T=0 Glauber tie
# rate (M2); 0.0 freezes ties, so at zero temperature a configuration with no
# energy-lowering flip is absorbing. Named presets: TIE_FLIP_PRESETS.
# Default 1.0 preserves the historical behaviour.
ISING_TIE_FLIP_DEFAULT: float = TIE_FLIP_PRESETS["standard"]

# ---------------------------------------------------------------------------
# Hamiltonian conventions (plan §3b)
# ---------------------------------------------------------------------------
# One symmetric coupling matrix J built once at setup; ΔE and the energy
# observable both derive from it. 'raw': J_ij = w_ij. 'sym': J_ij =
# w_ij/sqrt(deg_i deg_j) (normalized-adjacency Ising). 'avg': J_ij =
# w_ij (1/deg_i + 1/deg_j) (the per-site-averaged semantics today's
# calc_full_energy uses, done consistently).
ISING_COUPLING_NORMS: tuple[str, ...] = ("raw", "sym", "avg")
ISING_COUPLING_NORM_DEFAULT: str = "raw"

# Temperature (k_B = 1 units). 0.0 = zero-temperature (quench) dynamics —
# the historical constructor default.
ISING_T_DEFAULT: float = 0.0

# ---------------------------------------------------------------------------
# Observables (new scheme classes; ObservableSet names + file bases)
# ---------------------------------------------------------------------------
ISING_OBS_ENERGY: str = "energy"  # intensive per-site energy series
ISING_OBS_MAGN: str = "magn"  # per-spin magnetization series
ISING_OBS_SNAPSHOTS: str = "snapshots"  # (n_rec, N) spin trajectory

# Default recorded set — the legacy Metropolis sampler always recorded the
# energy + magnetization traces; snapshots are opt-in (I/O-heavy).
ISING_OBS_DEFAULT: tuple[str, ...] = (ISING_OBS_ENERGY, ISING_OBS_MAGN)

# On-disk file bases (``<base>_p=<pflip>[_suffix].bin`` under dynpath),
# matching the Voter/CP token conventions for the in-process backends. The C
# subprocess keeps its own frozen filenames (compat decision D-B5).
ISING_ENERGY_FBASE: str = "ene"
ISING_MAGN_FBASE: str = "m"
ISING_SNAPSHOTS_FBASE: str = "sout"

# Snapshot cadence + size guard (copied shape from ContactProcess): keep every
# k-th sweep when streaming, and refuse before writing if the subsampled
# trajectory would exceed the byte cap.
ISING_SNAPSHOT_EVERY_DEFAULT: int = 1
ISING_SNAPSHOT_MAX_BYTES: int = 512 * 1024 * 1024  # 512 MiB, as CP
