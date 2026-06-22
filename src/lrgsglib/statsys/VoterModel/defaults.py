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
#             majority. NOTE: this normalized power form is a project-chosen
#             parametrisation; it has no exact primary source among the refs
#             (citation PENDING) -- it is the nonlinear-voter family surveyed in
#             ref [4] Redner 2019, but no equation there matches it verbatim.
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
# balanced, i.e. lambda_min(signed Laplacian) <= FRUSTRATION_EIG_TOL
# (ref [3] Properties (ii)-(iii), p. 9). On a frustrated graph no absorbing
# state exists and the dynamics never freezes.
DEFAULT_ABSORBING_EVERY: int = 1     # check every N sweeps (1 = each sweep)
FRUSTRATION_EIG_TOL: float = 1e-10   # lambda_min(SL) above this => frustrated

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
# or b_ij = +1 ('rawspin'); one DSU engine serves both (the tracker only swaps
# the per-edge sign array it is handed). A single spin flip toggles the active
# state of EVERY incident edge, so the partition is maintained incrementally:
# union-find merges with the now-agreeing (formerly frustrated) neighbours in
# O(1), and a flip can only *split* the node's old domain when it had >= 2
# agreeing neighbours -- the only case that triggers a (bounded) local rescan.
CLUSTER_MODES: tuple[str, ...] = ("satisfied", "rawspin")
DEFAULT_CLUSTER_MODE: str = "satisfied"

# Integer codes shared with the native cluster tracker (LRGSG_clusters.c
# `rawspin` flag): 0 = satisfied (signed), 1 = rawspin (edge sign ignored).
CLUSTER_MODE_CODE: dict[str, int] = {"satisfied": 0, "rawspin": 1}

# On a potential split (>= 2 agreeing neighbours leave a domain) the tracker
# runs a bounded breadth-first scan from one departing neighbour: if it re-reaches
# all the others within CLUSTER_SPLIT_SCAN_CAP nodes the domain is intact (no
# split, O(cap)); otherwise it falls back to an exact recompute of just that one
# affected component (correctness guaranteed, O(component) only on a real split).
CLUSTER_SPLIT_SCAN_CAP: int = 256
