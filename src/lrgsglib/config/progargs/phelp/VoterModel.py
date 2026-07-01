"""Help strings for the L2D_VoterModel program + serializer CLIs.

Only the voter-specific strings live here; shared sweep / dispatch / slanzarv
help is reused from ``phelp.IsingDynamics`` and ``phelp.generic`` (imported by
``progargs.VoterModel``). The physics of each rule / schedule is documented on
the VoterModel methods and in ``statsys/VoterModel/defaults.py``.
"""

# --- Axis A: local update rule -------------------------------------------------
phelp_voter_rule = (
    "Local update rule: 'linear' (classic voter), 'majority' (local "
    "majority-vote), 'qvoter' (q-voter, sample q neighbours with replacement), "
    "or 'nonlinear' (power-alpha voter; alpha=1 reduces to linear)."
)
phelp_voter_q = (
    "q-voter only: number of neighbours sampled WITH replacement (>= 1)."
)
phelp_voter_eps = (
    "q-voter / nonlinear noise: flip probability when the sampled group is not "
    "unanimous (in [0, 1]); 0 keeps the rule deterministic-when-unanimous."
)
phelp_voter_alpha = (
    "nonlinear rule only: power exponent (>= 0). alpha=1 -> linear voter, "
    "alpha>1 reinforces the local majority, alpha<1 favours the minority."
)

# --- Axis B: update schedule ---------------------------------------------------
phelp_voter_upd_mode = (
    "Update schedule: 'asynchronous' (N single-node updates per sweep), "
    "'synchronous' (all nodes from a frozen snapshot), 'link' (edge-update), or "
    "'gillespie' (rejection-free CTMC). 'link' and 'gillespie' require "
    "rule='linear'."
)

# --- backend / init / time -----------------------------------------------------
phelp_voter_runlang = (
    "Execution backend: 'py' (reference), 'pb' (pybind in-process), 'np' "
    "(NumPy synchronous), 'cu' (CuPy GPU synchronous), or 'C'/'C0'/'C0S' "
    "(C subprocess). 'np'/'cu' force synchronous updates and cannot save the "
    "trajectory or track clusters."
)
phelp_voter_ic = (
    "Initial condition (e.g. 'uniform', 'random', 'ground_state_<idx>')."
)
phelp_voter_steps = (
    "Number of Monte Carlo sweeps (each sweep performs N single-node updates)."
)
phelp_voter_simref = (
    "Size-normalised time; steps = simref x N. Mutually exclusive with --steps."
)
phelp_voter_state_type = (
    "State encoding: 'bipolar' (-1/+1) or 'binary' (0/1)."
)
phelp_voter_freq = "Recording frequency for the history/observable sampling."
phelp_voter_num_log_samples = (
    "Log-spaced snapshot count for the C-subprocess backend (C0S mode)."
)
phelp_voter_seed = "RNG seed (auto-generated when omitted)."
phelp_voter_randstr = (
    "Append a random string to the run identifier so concurrent runs write to "
    "distinct output files."
)

# --- absorbing-state detection -------------------------------------------------
phelp_voter_absorbing_check = (
    "Detect and early-stop at a zero-frustration (absorbing) configuration. "
    "Omit to auto-select from the substrate."
)
phelp_voter_absorbing_every = "Check the absorbing condition every N sweeps."

# --- cluster-size distribution -------------------------------------------------
phelp_voter_cluster_mode = (
    "Active-edge predicate for the cluster-size distribution (when --save-cldist "
    "is set): 'satisfied' (signed domains, sign(w_ij)*s_i*s_j=+1) or 'rawspin' "
    "(s_i==s_j, edge sign ignored)."
)

# --- trajectory subsampling ----------------------------------------------------
phelp_voter_sout_every = (
    "Subsample the saved trajectory: keep every Nth recorded sweep (>= 1)."
)
phelp_voter_sout_nlog = (
    "Log-spaced snapshot count for the saved trajectory (overrides --sout-every)."
)
phelp_voter_sout_force_full = (
    "Force the full trajectory, overriding the 4 GiB soft cap fallback."
)

# --- observable save gates -----------------------------------------------------
phelp_voter_save_magn = "Record/persist the per-spin magnetization series (magn)."
phelp_voter_save_sout = "Record/persist the spin-configuration trajectory (sout)."
phelp_voter_save_cldist = (
    "Record/persist the cluster-size-distribution time series (cldist); enables "
    "cluster tracking."
)
phelp_voter_save_all = "Enable every observable (magn + sout + cldist)."

# --- serializer voter sweep axes ----------------------------------------------
phelp_voter_rule_list = (
    "Serializer sweep: one job per rule in this list (subset of "
    "{linear, majority, qvoter, nonlinear}). Omit to use the single --rule."
)
phelp_voter_upd_mode_list = (
    "Serializer sweep: one job per update schedule in this list. Omit to use the "
    "single --upd_mode."
)
phelp_voter_alpha_list = (
    "Serializer sweep: nonlinearity exponents (nonlinear rule). Omit for --alpha."
)
phelp_voter_alpha_linsp = (
    "Serializer sweep: semicolon-separated linspace tuples for alpha, e.g. "
    "'(0.5, 2.0, 7)'."
)
phelp_voter_eps_list = (
    "Serializer sweep: noise values eps (qvoter / nonlinear). Omit for --eps."
)
phelp_voter_eps_linsp = (
    "Serializer sweep: semicolon-separated linspace tuples for eps."
)
phelp_voter_q_list = (
    "Serializer sweep: q-voter neighbour counts. Omit to use the single --q."
)
