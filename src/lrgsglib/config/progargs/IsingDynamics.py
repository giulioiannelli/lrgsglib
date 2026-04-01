#
from .phelp.IsingDynamics import *
from .defs.IsingDynamics import *
#
from .Lattice2D import *
from .Lattice3D import *
from .ErdosRenyi import *
# general arguments
IsDyn_args: dict = {}
IsDyn_opt_args = {
    tuple(['-T', '--temperature']): {
        'help': phelp_T,
        'type': float,
        'default': 1.0,
        'dest': 'T',
    },
    tuple(['-fq', '--freq']): {
        'help': phelp_freq,
        'type': int,
        'default': DEFAULT_FREQ
    },
    tuple(['-ic', '--init_cond']): {
        'help': phelp_ic,
        'type': str,
        'default': DEFAULT_INIT_COND
    },
    tuple(['-nc', '--NoClust']): {
        'help': phelp_NoClust,
        'type': int,
        'default': DEFAULT_NOCLUST
    },
    tuple(['-os', '--out_suffix']): {
        'help': phelp_outsuffix,
        'type': str,
        'default': DEFAULT_OUTSFFX
    },
    tuple(['-rl', '--runlang']): {
        'help': phelp_runlang,
        'type': str,
        'default': DEFAULT_RUNLANG
    },
    tuple(['-ts', '--thrmsteps']): {
        'help': phelp_thrmsteps,
        'type': int,
        'default': DEFAULT_THRMSTEPS
    },
    tuple(['-es', '--eqstep']): {
        'help': phelp_eqstep,
        'type': int,
        'default': DEFAULT_EQSTEP
    },
    tuple(['-vl', '--val']): {
        'help': phelp_val,
        'type': str,
        'default': DEFAULT_VAL
    },
    tuple(['-rnds', '--randstr']): {
        'help': phelp_randstr_ising,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_randstr,
    },
    tuple(['--n_thermal']): {
        'help': phelp_n_thermal,
        'type': int,
        'default': DEFAULT_N_THERMAL,
    },
    tuple(['-sf', '--save_frequency']): {
        'help': phelp_save_frequency,
        'type': int,
        'default': DEFAULT_SAVE_FREQUENCY,
    },
    tuple(['--save-ene']): {
        'help': phelp_save_ene,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SAVE_ENE,
    },
    tuple(['--save-magn']): {
        'help': phelp_save_magn,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SAVE_MAGN,
    },
    tuple(['--save-sout']): {
        'help': phelp_save_sout,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SAVE_SOUT,
    },
    tuple(['--save-all']): {
        'help': phelp_save_all,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SAVE_ALL,
    },
}
# Simulated Annealing arguments
IsDyn_sa_args = {
    tuple(['--sa-mode']): {
        'help': phelp_sa_mode,
        'action': argparse.BooleanOptionalAction,
        'default': False,
    },
    tuple(['--T-init']): {
        'help': phelp_sa_T_init,
        'type': float,
        'default': DEFAULT_SA_T_INIT,
    },
    tuple(['--T-final']): {
        'help': phelp_sa_T_final,
        'type': float,
        'default': DEFAULT_SA_T_FINAL,
    },
    tuple(['--cooling-schedule']): {
        'help': phelp_sa_cooling_schedule,
        'type': str,
        'default': DEFAULT_SA_COOLING_SCHEDULE,
        'choices': ['linear', 'exponential', 'logarithmic'],
    },
    tuple(['--cooling-rate']): {
        'help': phelp_sa_cooling_rate,
        'type': float,
        'default': DEFAULT_SA_COOLING_RATE,
    },
    tuple(['--steps-per-T']): {
        'help': phelp_sa_steps_per_T,
        'type': int,
        'default': DEFAULT_SA_STEPS_PER_T,
    },
    tuple(['--n-temperatures']): {
        'help': phelp_sa_n_temperatures,
        'type': int,
        'default': DEFAULT_SA_N_TEMPERATURES,
    },
}
# Parallel Tempering arguments
IsDyn_pt_args = {
    tuple(['--pt-mode']): {
        'help': phelp_pt_mode,
        'action': argparse.BooleanOptionalAction,
        'default': False,
    },
    tuple(['--n-replicas']): {
        'help': phelp_pt_n_replicas,
        'type': int,
        'default': DEFAULT_PT_N_REPLICAS,
    },
    tuple(['--T-min']): {
        'help': phelp_pt_T_min,
        'type': float,
        'default': DEFAULT_PT_T_MIN,
    },
    tuple(['--T-max']): {
        'help': phelp_pt_T_max,
        'type': float,
        'default': DEFAULT_PT_T_MAX,
    },
    tuple(['--T-ladder-type']): {
        'help': phelp_pt_T_ladder_type,
        'type': str,
        'default': DEFAULT_PT_T_LADDER_TYPE,
        'choices': ['geometric', 'linear', 'custom'],
    },
    tuple(['--steps-per-exchange']): {
        'help': phelp_pt_steps_per_exchange,
        'type': int,
        'default': DEFAULT_PT_STEPS_PER_EXCHANGE,
    },
    tuple(['--n-exchanges']): {
        'help': phelp_pt_n_exchanges,
        'type': int,
        'default': DEFAULT_PT_N_EXCHANGES,
    },
}
# Topological algorithm arguments
IsDyn_topo_args = {
    tuple(['--topo-n-modes']): {
        'help': phelp_topo_n_modes,
        'type': int,
        'default': DEFAULT_TOPO_N_MODES,
    },
    tuple(['--topo-sigma-init']): {
        'help': phelp_topo_sigma_init,
        'type': float,
        'default': DEFAULT_TOPO_SIGMA_INIT,
    },
    tuple(['--topo-chunk-size']): {
        'help': phelp_topo_chunk_size,
        'type': int,
        'default': DEFAULT_TOPO_CHUNK_SIZE,
    },
    tuple(['--topo-polish']): {
        'help': phelp_topo_polish,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_TOPO_POLISH,
    },
    tuple(['--topo-polish-sweeps']): {
        'help': phelp_topo_polish_sweeps,
        'type': int,
        'default': DEFAULT_TOPO_POLISH_SWEEPS,
    },
    tuple(['--topo-tau']): {
        'help': phelp_topo_tau,
        'type': float,
        'default': DEFAULT_TOPO_TAU,
    },
    tuple(['--topo-field-strength']): {
        'help': phelp_topo_field_strength,
        'type': float,
        'default': DEFAULT_TOPO_FIELD_STRENGTH,
    },
}
# Cross-Entropy Method arguments
IsDyn_cem_args = {
    tuple(['--cem-iter']): {
        'help': phelp_cem_iter,
        'type': int,
        'default': DEFAULT_CEM_ITER,
    },
    tuple(['--cem-pop-size']): {
        'help': phelp_cem_pop_size,
        'type': int,
        'default': DEFAULT_CEM_POP_SIZE,
    },
    tuple(['--cem-elite-frac']): {
        'help': phelp_cem_elite_frac,
        'type': float,
        'default': DEFAULT_CEM_ELITE_FRAC,
    },
    tuple(['--cem-init-sigma']): {
        'help': phelp_cem_init_sigma,
        'type': float,
        'default': DEFAULT_CEM_INIT_SIGMA,
    },
    tuple(['--cem-smoothing']): {
        'help': phelp_cem_smoothing,
        'type': float,
        'default': DEFAULT_CEM_SMOOTHING,
    },
    tuple(['--cem-sigma-floor']): {
        'help': phelp_cem_sigma_floor,
        'type': float,
        'default': DEFAULT_CEM_SIGMA_FLOOR,
    },
    tuple(['--cem-sigma-ceiling']): {
        'help': phelp_cem_sigma_ceiling,
        'type': float,
        'default': DEFAULT_CEM_SIGMA_CEILING,
    },
    tuple(['--cem-restarts']): {
        'help': phelp_cem_restarts,
        'type': int,
        'default': DEFAULT_CEM_RESTARTS,
    },
    tuple(['--cem-greedy']): {
        'help': phelp_cem_greedy,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_CEM_GREEDY,
    },
    tuple(['--cem-greedy-sweeps']): {
        'help': phelp_cem_greedy_sweeps,
        'type': int,
        'default': DEFAULT_CEM_GREEDY_SWEEPS,
    },
}
IsDyn_srun_list_args = {
    tuple(['-s1', '--side1_list']): {
        'help': phelp_side1_list,
        'type': int,
        'nargs': '+',
        'default': DEFAULT_SIDE1_LIST
    },
    tuple(['-pFT', '--pflip_linsp']): {
        'help': phelp_pflip_linsp,
        'type': parse_multiple_linspace,
        'default': DEFAULT_PFLIP_LINSP
    },
    tuple(['-TT', '--Temp_linsp']): {
        'help': phelp_Temp_linsp,
        'type': parse_multiple_linspace,
        'default': DEFAULT_TEMP_LINSP
    },
}
IsDyn_srun_args = {
    tuple(['-nc', '--NoClust']): {
        'help': phelp_NoClust,
        'type': int,
        'default': DEFAULT_NOCLUST
    },
    tuple(['-mMB', '--slanzarv_minMB']): {
        'help': phelp_mMB,
        'type': int,
        'default': DEFAULT_mMB
    },
    tuple(['-MMB', '--slanzarv_maxMB']): {
        'help': phelp_MMB,
        'type': int,
        'default': DEFAULT_MMB
    },
    tuple(['--moretime']): {
        'help': phelp_moretime,
        'type': int,
        'default': DEFAULT_MORETIME
    },
    tuple(['--slanzarv_id']): {
        'help': phelp_slanzarv_id,
        'type': str,
        'default': DEFAULT_SLANZARV_ID
    },
}
# Lattice2D args
## names and descriptions
L2D_ISDYN_progname = 'L2D_IsingDynamics'
L2D_ISDYN_progname_shrt = 'L2DID'
L2D_ISDYN_description = f"""
    Computational resourses regarding the Ising Dynamics of 2D 
    lattices: {L2D_ISDYN_progname}.py
"""
L2D_ISDYN_srun_description = f"""Serialiser for {L2D_ISDYN_progname}.py"""
## arg parsers dict
L2D_ISDYN_args = {**L2D_args, **IsDyn_args}
L2D_ISDYN_action_args = {
    **action_args_dict,
    tuple(['--save-results']): {
        'help': phelp_save_results,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SAVE_RESULTS,
    },
}
L2D_ISDYN_opt_args = {**IsDyn_opt_args, **IsDyn_sa_args, **IsDyn_pt_args, **IsDyn_topo_args, **IsDyn_cem_args}
L2D_ISDYN_srun_opt_args = {**IsDyn_srun_list_args, **IsDyn_srun_args}
L2D_ISDYN_srun_action_args = {**srun_action_args}
# Lattice3D args
## names and descriptions
L3D_ISDYN_progname = 'L3D_IsingDynamics'
L3D_ISDYN_progname_shrt = 'L3DID'
L3D_ISDYN_description = f"""
    Computational resourses regarding the Ising Dynamics of 3D
    lattices: {L3D_ISDYN_progname}.py
"""
L3D_ISDYN_srun_description = f"""Serialiser for {L3D_ISDYN_progname}.py"""
## arg parsers dict
L3D_ISDYN_args = {**L3D_args, **IsDyn_args}
L3D_ISDYN_action_args = {
    **action_args_dict,
    tuple(['--save-results']): {
        'help': phelp_save_results,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SAVE_RESULTS,
    },
}
L3D_ISDYN_opt_args = {**IsDyn_opt_args, **IsDyn_sa_args, **IsDyn_pt_args, **IsDyn_topo_args, **IsDyn_cem_args}
L3D_ISDYN_srun_opt_args = {**IsDyn_srun_list_args, **IsDyn_srun_args}
L3D_ISDYN_srun_action_args = {**srun_action_args}
# ErdosRenyi args
## names and descriptions
ER_ISDYN_progname = 'ER_IsingDynamics'
ER_ISDYN_progname_shrt = 'ERID'
ER_ISDYN_description = f"""
    Computational resources regarding the Ising Dynamics on signed
    Erdos-Renyi graphs: {ER_ISDYN_progname}.py
"""
## arg parsers dict
ER_ISDYN_args = {**ER_args, **IsDyn_args}
ER_ISDYN_action_args = {**action_args_dict}
ER_ISDYN_opt_args = {**IsDyn_opt_args, **IsDyn_sa_args, **IsDyn_pt_args, **IsDyn_topo_args, **IsDyn_cem_args}