phelp_rrg_N = "Number of nodes N of the random regular graph"
phelp_rrg_k = "Degree k of every node (N*k even)"
phelp_rrg_seed = (
    "Root seed of the seed tree (the s= token of the run directory); "
    "another seed gives an independent replica"
)
phelp_rrg_sampler = (
    "Topology sampler: 'pairing' = exact uniform, 'native' = the engine's "
    "generator (NX native is biased), 'auto' = pairing for small k"
)
phelp_rrg_engine = "Graph engine that builds each topology: 'nx' or 'gt'"
phelp_tc_cell = (
    "Defect family: 'rand' = S (negative edges, control p = M-/M), "
    "'randXERR' = X (defect nodes with all edges negative, control q = "
    "node fraction)"
)
phelp_tc_ensemble = (
    "Disorder ensemble of S and X: 'canonical' = exact count "
    "floor(p*M + 1/2) (round half up), 'bernoulli' = independent draws"
)
phelp_tc_p_list = "Control values p (q for -c randXERR), space separated"
phelp_tc_p_linsp = (
    "Control values p as linspace tuples '(start, stop, num)' joined by ';'"
)
phelp_tc_sigma_list = "Coupling spreads sigma of -ew normal, space separated"
phelp_tc_sigma_linsp = (
    "Coupling spreads sigma as linspace tuples '(start, stop, num)' "
    "joined by ';'"
)
phelp_tc_rep_start = (
    "Index of the first rep (disorder realization); raise it to add "
    "statistics to an existing run"
)
phelp_tc_ctrl_chunk = "Control values per tile (0 = all)"
phelp_tc_rep_chunk = "Reps per tile (0 = all)"
phelp_tc_tile_id = (
    "Tile to run (0-based, rep-major); -1 runs every tile in this process"
)
phelp_tc_eig_tol = (
    "Lanczos tolerance (0 = machine precision); part of the run identity"
)
phelp_rrg_tag = (
    "Free tag of the run directory (the tag= token); it cannot contain _, "
    "= or /"
)
phelp_tc_gap_every = (
    "Reps with rep %% n == 0 also solve lambda_1 (gap, Davis-Kahan bound, "
    "e1); 0 = off"
)
# RRG_TransCluster_Serializer.py
phelp_rrg_N_list = "Numbers of nodes N; one run per (N, k)"
phelp_rrg_k_list = "Degrees k; one run per (N, k)"
phelp_tc_srun_actions = {
    "print": "Print the sbatch scripts; no side effects (the default)",
    "exec": (
        "Publish and check every cfg.json, then submit the arrays with "
        "sbatch; each one is logged to submissions/<job id>.json in its run "
        "directory"
    ),
    "status": "Count complete, partial and missing tiles per (N, k)",
    "collect": (
        "Aggregate the complete tiles into agg.npz and print P_inf, chi and "
        "the solve cost per control value"
    ),
}
phelp_tc_only_missing = (
    "Print or submit only the tiles without a complete file (running tiles "
    "count as missing: use it after the arrays have ended)"
)
phelp_tc_verify = (
    "Read every complete tile file instead of trusting its name: --status "
    "also counts unreadable and foreign tiles, and --only-missing "
    "resubmits the unreadable ones"
)
phelp_tc_array_concurrent = (
    "Tasks of one array running at once (the %%K). One exec submits an "
    "array per --array-max block and per (N, k), so up to their sum run at "
    "once"
)
phelp_tc_array_partition = "#SBATCH --partition (blank = cluster default)"
phelp_tc_array_time = "#SBATCH --time of each task (blank = cluster default)"
phelp_tc_array_max = (
    "Tasks per array; larger runs split into several arrays (keep it below "
    "the cluster's MaxArraySize)"
)
phelp_tc_mem_mb = (
    "#SBATCH --mem of each task in MB (0 = estimate from N and the engine)"
)
phelp_tc_boot = "Bootstrap resamples of the chi peak (0 = skip)"
phelp_tc_boot_seed = "Seed of the chi-peak bootstrap, in [0, 2**63)"
phelp_tc_peak_window = "Control window LO HI of the chi-peak fit (None = all)"
phelp_tc_legacy_txt = (
    "Also write the legacy 6-column ordParam.txt next to agg.npz; an "
    "existing one is refreshed by every --collect"
)
