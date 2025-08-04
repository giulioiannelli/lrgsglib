from .common import *
#
from .phelp.Topologies import *
from .phelp.Recon import *
#
from .defs.Topologies import *
from .defs.Recon import *
#
from .IsingDynamics import *
#
Recon_progName = 'Recon'
#
recon_opt_args_dict = {
    tuple(['-m', '--mode']): {
        'help': phelp_mode,
        'type': str,
        'default': DEFAULT_MODE
    },
    tuple(['-sFQ', '--save_frequency']):{
        'help': phelp_save_frequency,
        'type': int,
        'default': DEFAULT_SAVING_FREQUENCY,
    },
    tuple(['-mF', '--max_factor']): {
        'help': phelp_max_factor,
        'type': float,
        'default': DEFAULT_MAX_FACTOR
    },
    tuple(['-bS', '--basis_step']): {
        'help': phelp_basis_step,
        'type': int,
        'default': DEFAULT_BASIS_STEP
    },
    tuple(['-sib', '--start_idx_basis']): {
        'help': phelp_start_idx_basis,
        'type': int,
        'default': DEFAULT_START_IDX_BASIS
    },

}
# Lattice2D
## names and descriptions
L2D_Recon_progName = 'L2D_Recon'
L2D_Recon_progNameShrt = 'L2DRT'
L2D_Recon_description = f"""
    Computational resourses regarding the reconstruction of laplacian eigenmodes
    of 2D lattices: {L2D_Recon_progName}.py
"""
L2D_Recon_srun_description = f"""Serialiser for {L2D_Recon_progName}.py"""
## arg parsers dict
L2D_Recon_args = {**L2D_ISDYN_args}
L2D_Recon_optional_args_dict = {**recon_opt_args_dict}
L2D_Recon_action_args_dict = {**action_args_dict}
# Lattice3D
## names and descriptions
L3D_Recon_progName = 'L3D_Recon'
L3D_Recon_progNameShrt = 'L3DRT'
L3D_Recon_description = f"""
    Computational resourses regarding the reconstruction of laplacian eigenmodes
    of 3D lattices: {L3D_Recon_progName}.py
"""
L3D_Recon_srun_description = f"""Serialiser for {L3D_Recon_progName}.py"""
## arg parsers dict
L3D_Recon_args = {**L3D_ISDYN_args}
L3D_Recon_optional_args_dict = {**recon_opt_args_dict}
L3D_Recon_action_args_dict = {**action_args_dict}
#
L3D_Recon_srun_optional_args_dict = {
    tuple(['-dL', '--dim_list']): {
        'help': phelp_l3ddim_list,
        'type': int,
        'nargs': '+',
        'default': DEFAULT_L3DDIM_LIST
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
    **srun_opt_args
}