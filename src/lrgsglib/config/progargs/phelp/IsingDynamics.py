phelp_ic = "Initial condition for the Ising model"
phelp_runlang = (
    "Backend for running the Ising model. "
    "Options: py (Python Metropolis), pb_met (pybind11 Metropolis), "
    "pb_sa (pybind11 SA), pb_pt (pybind11 PT), "
    "cu_met/cu_sa/cu_pt (CuPy GPU), "
    "wolff/sw (cluster algorithms), "
    "C<digit><letters> (C subprocess, NX only): "
    "digit: 0=Metropolis, 1=SA, 2=PT; "
    "letters: E=energy/magn, S=snapshots, K=cluster, V=eigvec, H=hfield"
)
phelp_thrmsteps = "Number of thermalization steps"
phelp_eqstep = "Number of equilibrium/optimization MC steps (1 step = N proposals)"
phelp_randstr = "Random string for the output files"
phelp_quench_id = (
    "Run exactly one quench with this 1-indexed id (for SLURM array "
    "dispatch). If -1 (default), loop 1..number_of_averages internally."
)
# Array-dispatch help text (serializer only)
phelp_dispatch = (
    "Dispatch mode: 'slanzarv' (default, 1 job per parameter point) or "
    "'array' (SLURM --array with 1 task per quench realization)."
)
phelp_topo_n_modes_list = (
    "Sweep topo_n_modes values; one submission per M. Overrides "
    "--topo-n-modes when provided."
)
phelp_cem_iter_list = (
    "Per-M cem_iter values; length must equal -mL length. Overrides "
    "--cem-iter when provided."
)
phelp_array_concurrent = (
    "SLURM --array concurrency throttle (the %%N in 'array=1-N%%M')."
)
phelp_array_partition = (
    "#SBATCH --partition= value for array jobs. Blank uses cluster default."
)
phelp_array_time = (
    "#SBATCH --time= value for array jobs (e.g. 08:00:00). Blank = default."
)
# Simulated Annealing help text
phelp_sa_mode = "Enable simulated annealing mode"
phelp_sa_T_init = "Initial temperature for SA"
phelp_sa_T_final = "Final temperature for SA"
phelp_sa_cooling_schedule = "Temperature cooling schedule"
phelp_sa_cooling_rate = "Cooling rate (for exponential schedule)"
phelp_sa_steps_per_T = "MC steps at each temperature"
phelp_sa_n_temperatures = "Number of temperature steps"
# Parallel Tempering help text
phelp_pt_mode = "Enable parallel tempering mode"
phelp_pt_n_replicas = "Number of temperature replicas"
phelp_pt_T_min = "Minimum temperature in ladder"
phelp_pt_T_max = "Maximum temperature in ladder"
phelp_pt_T_ladder_type = "Type of temperature ladder"
phelp_pt_steps_per_exchange = "MC sweeps between exchange attempts"
phelp_pt_n_exchanges = "Total number of exchange rounds"
# Topological algorithm help text
phelp_topo_n_modes = "Number of Laplacian eigenvectors for topological algorithms"
phelp_topo_sigma_init = "Initial proposal width for topo_met"
phelp_topo_chunk_size = "Steps between adaptive sigma adjustments"
phelp_topo_polish = "Enable T=0 greedy polish after each chunk"
phelp_topo_polish_sweeps = "Max sweeps per greedy polish"
phelp_topo_tau = "Softmax temperature for TFCA field weights"
phelp_topo_field_strength = "Overall field magnitude scale for TFCA"
# Cross-Entropy Method help text
phelp_cem_iter = "Number of CEM iterations per restart"
phelp_cem_pop_size = "Population size (K) for CEM sampling"
phelp_cem_elite_frac = "Fraction of population used as elite set"
phelp_cem_init_sigma = "Initial std deviation for CEM coefficient sampling"
phelp_cem_smoothing = "Exponential smoothing factor for CEM distribution updates"
phelp_cem_sigma_floor = "Minimum allowed sigma (prevents premature convergence)"
phelp_cem_sigma_ceiling = "Maximum allowed sigma (prevents divergence)"
phelp_cem_restarts = "Number of independent CEM restarts (best-of-R)"
phelp_cem_greedy = "Apply T=0 greedy quench to each CEM sample"
phelp_cem_greedy_sweeps = "Max sweeps per greedy quench within CEM"
# Thermal averaging help text
phelp_n_thermal = "Number of independent thermal runs per disorder realization"
# Result persistence help text
phelp_save_results = "Master switch: enable saving results to NPZ files"
phelp_save_frequency = "Checkpoint every N thermal runs (0 = save only at end)"
phelp_save_ene = "Save energy trajectories"
phelp_save_magn = "Save magnetization trajectories"
phelp_save_sout = "Save spin snapshots at log-spaced time points"
phelp_save_all = "Save all observables (ene + magn + sout + topo coeffs)"
