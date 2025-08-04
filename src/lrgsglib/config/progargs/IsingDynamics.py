#
from .phelp.IsingDynamics import *
from .defs.IsingDynamics import *
#
from .Lattice2D import *
from .Lattice3D import *
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
    tuple(['-rnds', '--rnd_str']): {
        'help': phelp_rnd_str,
        'type': bool,
        'default': DEFAULT_RND_STR
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
L2D_ISDYN_opt_args = {**IsDyn_opt_args}
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
L3D_ISDYN_opt_args = {**IsDyn_opt_args}
L3D_ISDYN_srun_opt_args = {**IsDyn_srun_list_args, **IsDyn_srun_args}
L3D_ISDYN_srun_action_args = {**srun_action_args}