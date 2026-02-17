phelp_ic = "Initial condition for the Ising model"
phelp_runlang = (
    "Backend for running the Ising model. "
    "Options: py (Python Metropolis), pb_met (pybind11 Metropolis), "
    "pb_sa (pybind11 SA), pb_pt (pybind11 PT), "
    "cu_met/cu_sa/cu_pt (CuPy GPU), "
    "wolff/sw (cluster algorithms), "
    "C1/C1b/C3/C3b/C4/C4b (C subprocess, NX only)"
)
phelp_thrmsteps = "Number of thermalization steps"
phelp_randstr = "Random string for the output files"
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
# Thermal averaging help text
phelp_n_thermal = "Number of independent thermal runs per disorder realization"
