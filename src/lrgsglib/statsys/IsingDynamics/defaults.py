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
# calc_full_energy uses, done consistently). The vocabulary is single-sourced
# in the shared coupling builder (statsys._csr).
ISING_COUPLING_NORMS: tuple[str, ...] = COUPLING_NORMS
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

# ---------------------------------------------------------------------------
# Simulated annealing (IsingSimulatedAnnealing; legacy ctor values preserved)
# ---------------------------------------------------------------------------
# Cooling schedules: 'linear' interpolates T_init → T_final over the stages;
# 'exponential' is T_init · rate^k (T_final unused, as in the legacy code);
# 'logarithmic' is T_init / ln(2 + k).
ISING_SA_SCHEDULES: tuple[str, ...] = (
    "linear",
    "exponential",
    "logarithmic",
)
ISING_SA_SCHEDULE_DEFAULT: str = "exponential"
ISING_SA_T_INIT_DEFAULT: float = 10.0
ISING_SA_T_FINAL_DEFAULT: float = 0.01
ISING_SA_COOLING_RATE_DEFAULT: float = 0.95
ISING_SA_N_TEMPERATURES_DEFAULT: int = 100
ISING_SA_STEPS_PER_T_DEFAULT: int = 100

# Legacy TFCA ("topological field-cooling annealing") reproduction: the
# legacy algorithm is a sudden quench at this constant temperature for
# n_temperatures * steps_per_T sweeps under the spectral field; the strangler
# facade maps it to IsingSimulatedAnnealing(field='spectral',
# T_schedule=full(n_temperatures, ISING_TFCA_QUENCH_T)).
ISING_TFCA_QUENCH_T: float = 1e-3

# External-field construction axis on the annealing scheme (plan §10):
# 'uniform' uses whatever field vector the user passed (default zero);
# 'spectral' builds the softmax-weighted eigenvector field at init (the
# legacy TFCA field) via statsys._spectral.build_spectral_field.
ISING_FIELD_MODES: tuple[str, ...] = ("uniform", "spectral")
ISING_FIELD_MODE_DEFAULT: str = "uniform"

# ---------------------------------------------------------------------------
# Spectral capability (move='spectral', field='spectral', IsingCEM subspace)
# ---------------------------------------------------------------------------
ISING_SPECTRAL_N_MODES_DEFAULT: int = 40
# Softmax temperature + magnitude of the spectral (TFCA) field.
ISING_SPECTRAL_TAU_DEFAULT: float = 1.0
ISING_SPECTRAL_FIELD_STRENGTH_DEFAULT: float = 1.0
# Spectral Metropolis move: initial per-mode proposal width, adaptation
# cadence (proposals between σ updates), and the optional T=0 greedy polish.
ISING_SPECTRAL_SIGMA_INIT_DEFAULT: float = 0.15
ISING_SPECTRAL_CHUNK_SIZE_DEFAULT: int = 5000
ISING_SPECTRAL_POLISH_DEFAULT: bool = True
ISING_SPECTRAL_POLISH_SWEEPS_DEFAULT: int = 50

# ---------------------------------------------------------------------------
# Parallel tempering (IsingParallelTempering; legacy ctor values preserved)
# ---------------------------------------------------------------------------
ISING_PT_N_REPLICAS_DEFAULT: int = 8
ISING_PT_T_MIN_DEFAULT: float = 0.5
ISING_PT_T_MAX_DEFAULT: float = 5.0
# Ladder spacings: 'geometric' = T_min (T_max/T_min)^(r/(R−1)) (uniform in
# ln T — the standard choice); 'linear' = uniform in T.
ISING_PT_LADDERS: tuple[str, ...] = ("geometric", "linear")
ISING_PT_LADDER_DEFAULT: str = "geometric"
ISING_PT_STEPS_PER_EXCHANGE_DEFAULT: int = 10
ISING_PT_N_EXCHANGES_DEFAULT: int = 1000

# ---------------------------------------------------------------------------
# Cross-entropy method (IsingCEM; legacy ctor values preserved)
# ---------------------------------------------------------------------------
ISING_CEM_ITER_DEFAULT: int = 30
ISING_CEM_POP_SIZE_DEFAULT: int = 128
ISING_CEM_ELITE_FRAC_DEFAULT: float = 0.2
ISING_CEM_INIT_SIGMA_DEFAULT: float = 1.2
ISING_CEM_SMOOTHING_DEFAULT: float = 0.6
ISING_CEM_SIGMA_FLOOR_DEFAULT: float = 1e-3
ISING_CEM_SIGMA_CEILING_DEFAULT: float = 5.0
ISING_CEM_RESTARTS_DEFAULT: int = 10
ISING_CEM_GREEDY_DEFAULT: bool = True
ISING_CEM_GREEDY_SWEEPS_DEFAULT: int = 120
