from .common import *
from .phelp.SCSNN import *
from .defs.SCSNN import *

SCSGeneralized_args = {
    'N': {
        'help': phelp_scs_N,
        'type': int,
    },
}

SCSGeneralized_optional_args_dict = {
    tuple(['--gamma']): {
        'help': phelp_scs_gamma,
        'type': float,
        'default': DEFAULT_SCS_NN_GAMMA,
    },
    tuple(['--J']): {
        'help': phelp_scs_J,
        'type': float,
        'default': DEFAULT_SCS_NN_J,
    },
    tuple(['--g']): {
        'help': phelp_scs_g,
        'type': float,
        'default': DEFAULT_SCS_NN_G,
    },
    tuple(['--diagonal']): {
        'help': phelp_scs_diagonal,
        'type': str,
        'default': DEFAULT_SCS_NN_DIAGONAL,
    },
    tuple(['-wd', '--workdir']): {
        'help': phelp_scs_workdir,
        'type': str,
        'default': DEFAULT_SCS_NN_WORKDIR,
    },
    tuple(['--seed']): {
        'help': phelp_scs_seed,
        'type': int,
        'default': None,
    },
}

SCSGeneralized_action_args_dict = {**action_args_dict}
