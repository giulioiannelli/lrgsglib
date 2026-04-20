"""argparse help strings for the SignedRW walker CLI."""

phelp_srw_rule = (
    "Walker rule: 'absorb' (D3 — dies on killing edge), 'kill' (D2 — dies "
    "with prob 1/n_kill at current node), 'sticky' (D1 — P(move)=2^(-n_cross))."
)
phelp_srw_mode = (
    "Program mode: 'single' runs one walker ensemble and saves per-walker "
    "observables; 'pair' runs two independent ensembles with distinct start "
    "configs and additionally saves their overlap metrics."
)
phelp_srw_n_walkers = "Number of independent walkers in an ensemble."
phelp_srw_coverage = (
    "Coverage threshold — each walker terminates once its unique-visit "
    "fraction reaches this value."
)
phelp_srw_x_node = (
    "Unfrustrated-X-node behaviour: 'reflect' (library default — "
    "unfrustrated negative edges reflect the walker) or 'absorb' (any "
    "negative-edge crossing kills, matches the FrustratedRW.ipynb notebook)."
)
phelp_srw_max_n_cross = "Sticky-rule only — death threshold on frustrated crossings."
phelp_srw_start_a = (
    "Walker-A start protocol: 'random' (uniform over nodes), 'fixed' "
    "(uses --start-a-node), 'center' (lattice-aware, (s//2)*(s+1) index)."
)
phelp_srw_start_b = "Walker-B start protocol (pair mode only; same choices as --start-a)."
phelp_srw_start_node_a = "Node index for walker A when --start-a=fixed."
phelp_srw_start_node_b = "Node index for walker B when --start-b=fixed."
phelp_srw_seed_a = "RNG seed for walker A."
phelp_srw_seed_b = "RNG seed for walker B (pair mode only)."
phelp_srw_store_trajectory = "Save full per-step trajectories (memory-heavy)."
phelp_srw_store_per_walker_visits = (
    "Save per-walker (n_walkers, N) visit-count arrays (memory-heavy)."
)
phelp_srw_runlang = "Execution backend ('py' is the only option in Phase 1)."
