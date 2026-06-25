"""Defaults and enumerations for the VoterModel rule family and update schedules.

Centralises every tunable for the Phase-3 voter rule family and sampler axis so
no magic numbers leak into the dynamics code (project hard-ban #1). The physics
of each rule is documented (with primary-source citations) on the VoterModel
methods that consume these constants.

All cited models are defined on **unsigned** graphs in their primary sources;
the library's signed-substrate generalisation -- the effective neighbour
opinion ``sigma_j = sign(w_ij) * s_j`` (a negative edge gives an anti-copy) --
is a project extension, **not** attributable to those papers.

References
----------
.. [1] Castellano, Fortunato, Loreto, Rev. Mod. Phys. 81, 591 (2009).
       voter_dynamics/castellano2009social.pdf -- linear voter (Eq. 5, p. 599),
       local majority-vote (Sec. III.B.1, p. 600), Galam majority rule
       (Sec. III.C, p. 602), link-update (Sec. III.B.3, p. 601).
.. [2] Castellano, Munoz, Pastor-Satorras, Phys. Rev. E 80, 041129 (2009).
       voter_dynamics/castellano2009nonlinear.pdf -- non-linear q-voter,
       flip rate f(x,q)=x^q+eps[1-x^q-(1-x)^q] (Eq. 1, Sec. II, p. 2).
.. [3] Iannelli, Villegas, Gili, Gabrielli, arXiv:2504.00144 (2026).
       iannelli2025topological.pdf -- signed Laplacian L_bar = D_bar - A
       (Eq. 4, p. 9); lambda_min(L_bar)=0 iff balanced, >0 iff frustrated,
       and lambda_min measures frustration (Properties (i)-(iii), p. 9).
.. [4] Redner, C. R. Physique 20, 275 (2019). voter_dynamics/redner2019reality.pdf
       -- mini-review of reality-inspired voter extensions (survey context for the
       q-voter / nonlinear-voter / majority-vote families).
.. [5] Yang, Wang, Lai, Wang, Phys. Lett. A 376, 282 (2012); arXiv:1008.0901.
       voter_dynamics/yang2012nonlinear.pdf -- nonlinear (power-alpha) voter:
       agent i selects opinion +1 w.p. p_+ = n_+^alpha / (n_+^alpha + n_-^alpha)
       (Eq. 1, Sec. II, p. 1), with n_+/n_- the counts among agent i AND its
       nearest neighbours; alpha=1 -> voter, alpha->inf -> local majority, alpha<1
       -> minority-favouring.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Axis A -- local update rule (``rule=``)
# ---------------------------------------------------------------------------
# linear    : copy sign(w_ij)*s_j of one uniformly chosen neighbour (classic
#             voter, ref [1] Eq. 5).
# majority  : adopt sign(sum_j w_ij s_j) with a random tie-break (local
#             majority-vote, ref [1] Sec. III.B.1 p. 600).
# qvoter    : sample ``q`` neighbours WITH replacement; if their signed opinions
#             are unanimous the node copies them, otherwise it flips with
#             probability ``eps`` (non-linear q-voter, ref [2] Sec. II; sampling
#             with repetition is explicit in the paper).
# nonlinear : reinforcement voter -- adopt opinion sigma with probability
#             P(sigma) = f_sigma^alpha / (f_+^alpha + f_-^alpha), where f_sigma
#             is the fraction of neighbours holding the signed opinion sigma.
#             alpha=1 recovers the linear voter; alpha>1 reinforces the local
#             majority; alpha<1 favours the minority. This is the nonlinear
#             power-alpha voter of ref [5] Yang et al. 2012 (Eq. 1:
#             p_+ = n_+^alpha / (n_+^alpha + n_-^alpha); alpha=1 -> voter,
#             alpha->inf -> local majority). CAVEAT: Yang et al. count n_+/n_-
#             over agent i AND its neighbours (focal site INCLUDED); this repo
#             normalises over neighbour fractions f_sigma (focal site EXCLUDED).
#             The power form is identical (the degree normalisation cancels in the
#             ratio); the focal-site inclusion is the only difference from [5].
#             ``eps`` does NOT enter this rule (it is a q-voter parameter).
VOTER_RULES: tuple[str, ...] = ("linear", "majority", "qvoter", "nonlinear")
DEFAULT_RULE: str = "linear"

# Per-rule parameter defaults. Chosen identity-preserving where it matters:
# ``alpha=1`` makes ``nonlinear`` coincide with the linear voter, and the noise
# ``eps=0`` makes the rules deterministic-when-unanimous.
DEFAULT_QVOTER_Q: int = 2           # neighbours sampled with replacement (q-voter)
DEFAULT_NOISE_EPS: float = 0.0      # flip prob. when not unanimous (qvoter / nonlinear)
DEFAULT_NONLIN_ALPHA: float = 1.0   # nonlinearity exponent (alpha=1 -> linear voter)

# ---------------------------------------------------------------------------
# Axis B -- update schedule / sampler (``upd_mode=``)
# ---------------------------------------------------------------------------
# asynchronous : N single-node updates per sweep, nodes drawn WITH replacement
#                (random-sequential; matches the C kernel and the q-voter
#                t -> t + 1/N convention).
# synchronous  : every node updated simultaneously from a frozen snapshot of the
#                previous sweep (double-buffer).
# link         : edge-update (ref [1] Sec. III.B.3 p. 601, Suchecki 2005a) -- a
#                random edge is picked and one randomly chosen endpoint copies
#                the other; intrinsically a copy, so defined for rule='linear'.
# gillespie    : rejection-free continuous-time MC (BKL/kinetic MC). Each node
#                attempts at rate 1 (so t=1 is one sweep = N attempts) and, given
#                an attempt, copies a random signed neighbour -- so its node-flip
#                rate is r_i = f_i / deg_i with f_i the number of frustrated edges
#                incident to i. Flips are sampled directly (node i w.p. r_i/R,
#                dt ~ Exp(R), R = sum_i r_i), skipping null events; the run freezes
#                exactly when R=0 (no frustrated edge = absorbing). Statistically
#                equivalent to the asynchronous linear voter, with a large speedup
#                near consensus. Intrinsically a copy rule, so rule='linear' only.
#                Implemented on all backends (Python, C subprocess, pybind) via
#                the shared _ccore CTMC kernel (LRGSG_ctmc.{c,h}).
VOTER_UPD_MODES: tuple[str, ...] = (
    "asynchronous", "synchronous", "link", "gillespie",
)
VOTER_UPD_MODES_PLANNED: tuple[str, ...] = ()
DEFAULT_UPD_MODE: str = "asynchronous"

# Update schedules that only make sense for a copy rule.
LINK_ONLY_RULES: frozenset[str] = frozenset({"linear"})
GILLESPIE_RULES: frozenset[str] = frozenset({"linear"})

# Update schedules implemented in the Python backend but NOT yet in the native
# (C subprocess / pybind) backends -- those raise rather than silently lie.
# (Empty now that gillespie has a native CTMC kernel; kept for future modes.)
PY_ONLY_UPD_MODES: frozenset[str] = frozenset()

# ---------------------------------------------------------------------------
# Integer codes shared with the C / pybind backends.
# MUST stay in sync with the voter_rule_t / voter_upd_t enums in
# ccore/LRGSG_vm.h.
# ---------------------------------------------------------------------------
RULE_CODE: dict[str, int] = {
    "linear": 0,
    "majority": 1,
    "qvoter": 2,
    "nonlinear": 3,
}
UPD_MODE_CODE: dict[str, int] = {
    "asynchronous": 0,
    "synchronous": 1,
    "link": 2,
    "gillespie": 3,
}

# ---------------------------------------------------------------------------
# Absorbing-state detection (signed substrate)
# ---------------------------------------------------------------------------
# A voter configuration is absorbing iff it has zero frustrated edges
# (w_ij s_i s_j > 0 for every edge) -- consensus on an unsigned graph, the
# signed ground state on a balanced graph. Such a state exists iff the graph is
# balanced (lambda_min(signed Laplacian) = 0; ref [3] Properties (ii)-(iii),
# p. 9). Detection here is the exact integer frustrated-edge count
# (VoterModel.count_frustrated_edges); on a frustrated graph no absorbing state
# exists and the dynamics never freezes.
DEFAULT_ABSORBING_EVERY: int = 1     # check every N sweeps (1 = each sweep)

# Default order parameter for the ensemble susceptibility helper
# (VoterModel.order_parameter_susceptibility); one of the per-run scalar
# observables {"largest_domain_fraction", "interface_density"}.
VOTER_DEFAULT_OP: str = "largest_domain_fraction"

# ---------------------------------------------------------------------------
# Cluster-size distribution observable (domains of aligned signed opinion)
# ---------------------------------------------------------------------------
# A *cluster* is a connected component of the subgraph of *active* edges. The
# activation predicate depends on ``cluster_mode``:
#   satisfied : edge (i,j) active iff sign(w_ij) * s_i * s_j = +1 -- a domain on
#               the SIGNED substrate (a negative edge "agrees" when the spins are
#               opposite). Reduces to same-opinion domains on an unsigned graph.
#   rawspin   : edge (i,j) active iff s_i == s_j (edge sign ignored) -- domains of
#               equal raw spin value.
# Both predicates are b_ij * s_i * s_j > 0 with b_ij = sign(w_ij) ('satisfied')
# or b_ij = +1 ('rawspin'); one connected-components pass serves both (only the
# per-edge sign array changes). The distribution is recomputed from scratch at
# each recorded sweep (O(N + E)); with ~N single-spin events per recorded sweep
# this matches the cost of incremental maintenance while staying trivially
# correct (it IS the connected-components definition).
CLUSTER_MODES: tuple[str, ...] = ("satisfied", "rawspin")
DEFAULT_CLUSTER_MODE: str = "satisfied"

# Integer codes shared with the native cluster recompute (LRGSG_clusters.c
# `rawspin` flag): 0 = satisfied (signed), 1 = rawspin (edge sign ignored).
CLUSTER_MODE_CODE: dict[str, int] = {"satisfied": 0, "rawspin": 1}

# ---------------------------------------------------------------------------
# On-disk persistence of observables (``savedisk``; see VoterModel._persist_*)
# ---------------------------------------------------------------------------
# When ``savedisk=True`` (default) an *enabled* observable is written under the
# run's ``dynpath`` (``data/<lattice>/voter/N=<n>/``) and the in-RAM attribute
# becomes a lazy, disk-backed view (see VoterModel.magn / s_t / cluster_dist).
# File bases passed to ``sg.get_p_fname`` (filename ``<base>_p=<pflip>[_suf].<ext>``,
# mirroring the C-backend convention so an analysis program can stream big output
# to disk and load it back only when needed):
VOTER_MAGN_FBASE: str = "m"         # per-spin magnetization series (float64 .bin)
VOTER_SOUT_FBASE: str = "sout"      # spin-configuration trajectory (int8 .bin, (n_rec, N))
VOTER_CLDIST_FBASE: str = "cldist"  # cluster-size-distribution time series (.npz)
# Logical names of the three observables in VoterModel.observables (the keys for
# ``observables[...]`` and ``output_sizes()``; the magn name is "magn", NOT the
# "m" file base above).
VOTER_OBS_MAGN: str = "magn"
VOTER_OBS_SOUT: str = "sout"
VOTER_OBS_CLDIST: str = "cldist"
# NOTE: ``cldist`` is DISTINCT from Ising's ``outcl{i}`` (time-averaged
# magnetization moments of pre-defined eigenvector clusters); the voter artifact
# is a per-sweep histogram {component_size: count} of emergent active-edge
# connected components -- a time series of domain-size distributions.

# ``sout`` (the spin trajectory) is the only observable whose footprint scales
# with steps*N, so it is the one that can spill terabytes on extensive runs. It
# is subsampled to every ``sout_every`` recorded sweeps; if the precomputed size
# (n_rec * N int8 bytes) would exceed ``VOTER_SOUT_MAX_BYTES`` the run does NOT
# raise -- it falls back to ``VOTER_SOUT_NLOG_DEFAULT`` LOG-SPACED snapshots (the
# in-process analogue of the C ``nSampleLog``). Override with ``sout_nlog=K``
# (explicit count) or ``sout_force_full=True`` (force the full trajectory). The
# C subprocess streams its own snapshots, sampled by ``nSampleLog``.
DEFAULT_SOUT_EVERY: int = 1                       # 1 = every recorded sweep (full trajectory)
VOTER_SOUT_MAX_BYTES: int = 4096 * 1024 * 1024    # 4 GiB soft cap for the sout trajectory
VOTER_SOUT_NLOG_DEFAULT: int = 1000               # log-spaced snapshots used on fallback

# Registry key under which VoterModel's solver backends are registered in
# ``statsys._solver_engine`` (single-sourced so the model and its solvers agree).
VOTER_SOLVER_NAME: str = "voter"
