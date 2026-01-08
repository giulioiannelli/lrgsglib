from .common import *

from .phelp.SlaplSpect import *
from .defs.SlaplSpect import *

from .Lattice2D import *
from .MultiplicativeCascade import *

#
# Common SlaplSpect parameters (apply to all graph types)
#
SlaplSpect_progName = 'SlaplSpect'

# Generic spectral analysis parameters (shared across L2D, MCG, etc.)
SlaplSpect_common_optional_args_dict = {
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
    tuple(['-sf', '--save_frequency']): {
        'help': phelp_data_save_freq,
        'type': int,
        'default': DEFAULT_DATA_SAVE_FREQ
    },
    tuple(['--backend']): {
        'help': phelp_slaplspect_backend,
        'type': str,
        'default': DEFAULT_SLAPLSPECT_BACKEND,
        'choices': ['numpy', 'scipy', 'cupy']
    },
    tuple(['--keep-sparse']): {
        'help': phelp_keep_sparse,
        'type': lambda x: None if x.lower() == 'none' else x.lower() == 'true',
        'default': DEFAULT_KEEP_SPARSE,
        'metavar': 'BOOL|none'
    }
}

#
# L2D_SlaplSpect (Lattice2D spectral analysis)
#
L2D_SlaplSpect_progName = 'L2D_SlaplSpect'
L2D_SlaplSpect_progNameShrt = 'L2DSS'
L2D_SlaplSpect_description = f"""
    Computational resourses regarding the Signed Laplacian spectrum of 2D
    lattices: {L2D_SlaplSpect_progName}.py
"""

L2D_SlaplSpect_args = {**L2D_args}
L2D_SlaplSpect_optional_args_dict = {
    **SlaplSpect_common_optional_args_dict,  # Inherit common spectral parameters
    tuple(['-m', '--mode']): {
        'help': phelp_l2dsspect_mode,
        'type': str,
        'default': DEFAULT_L2DSSPECT_MODE
    }
}
L2D_SlaplSpect_action_args_dict = {}

#
# MCG_SlaplSpect (MultiplicativeCascade spectral analysis)
#
MCG_SlaplSpect_progName = 'MCG_SlaplSpect'
MCG_SlaplSpect_progNameShrt = 'MCGSS'
MCG_SlaplSpect_description = f"""
    Computational resources regarding the Signed Laplacian spectrum of
    MultiplicativeCascadeGraph: {MCG_SlaplSpect_progName}.py
"""

MCG_SlaplSpect_args = {**MC_args}  # Inherits MC graph parameters
MCG_SlaplSpect_optional_args_dict = {
    **SlaplSpect_common_optional_args_dict,  # Inherit common spectral parameters
    tuple(['-m', '--mode']): {
        'help': phelp_mcgsspect_mode,
        'type': str,
        'default': DEFAULT_MCGSSPECT_MODE
    }
}
MCG_SlaplSpect_action_args_dict = {}
