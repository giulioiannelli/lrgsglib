#
from .phelp.IsingDynamics import *
from .defs.IsingDynamics import *
#
from .Lattice2D import *
from .Lattice3D import *
from .ErdosRenyi import *
# general arguments
IsDyn_args = {
    'T': {
        'help': phelp_T,
        'type': float
    }
}
IsDyn_opt_args = {
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
L2D_ISDYN_action_args = {**action_args_dict}
L2D_ISDYN_opt_args = {**IsDyn_opt_args, **IsDyn_sa_args, **IsDyn_pt_args}
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
L3D_ISDYN_action_args = {**action_args_dict}
L3D_ISDYN_opt_args = {**IsDyn_opt_args, **IsDyn_sa_args, **IsDyn_pt_args}
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
ER_ISDYN_opt_args = {**IsDyn_opt_args, **IsDyn_sa_args, **IsDyn_pt_args}