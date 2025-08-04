from .common import *

from .phelp.SlaplSpect import *
from .defs.SlaplSpect import *

from .Lattice2D import *
#
SlaplSpect_progName = 'SlaplSpect'
## names and descriptions
L2D_SlaplSpect_progName = 'L2D_SlaplSpect'
L2D_SlaplSpect_progNameShrt = 'L2DSS'
L2D_SlaplSpect_description = f"""
    Computational resourses regarding the Signed Laplacian spectrum of 2D 
    lattices: {L2D_SlaplSpect_progName}.py
"""
## arg parsers dict
L2D_SlaplSpect_args = {**L2D_args}
L2D_SlaplSpect_optional_args_dict = {
    tuple(['-bc', '--bins_count']): {
        'help': phelp_binsc,
        'type': int,
        'default': DEFAULT_BINSC
    },
    tuple(['-em', '--eigen_mode']): {
        'help': phelp_eigMode,
        'type': str,
        'default': DEFAULT_EIGMODE
    },
    tuple(['-hm', '--howmany']): {
        'help': phelp_howmany_eigs,
        'type': int,
        'default': DEFAULT_HOWMANY_EIGS
    },
    tuple(['-m', '--mode']): {
        'help': phelp_l2dsspect_mode,
        'type': str,
        'default': DEFAULT_L2DSSPECT_MODE
    },
    tuple(['-prd', '--period']): {
        'help': phelp_data_save_freq,
        'type': int,
        'default': DEFAULT_DATA_SAVE_FREQ
    }
}
L2D_SlaplSpect_action_args_dict = {}