from .common import *
from .phelp.SCSGeneralized import *
from .defs.SCSGeneralized import *

SCSGeneralized_args = {
    'N': {
        'help': phelp_scs_N,
        'type': int,
    },
    'gamma': {
        'help': phelp_scs_gamma,
        'type': float,
    },
}

SCSGeneralized_optional_args_dict = {
    tuple(['--J']): {
        'help': phelp_scs_J,
        'type': float,
        'default': DEFAULT_SCS_GENERALIZED_J,
    },
    tuple(['--g']): {
        'help': phelp_scs_g,
        'type': float,
        'default': DEFAULT_SCS_GENERALIZED_G,
    },
    tuple(['--diagonal']): {
        'help': phelp_scs_diagonal,
        'type': str,
        'default': DEFAULT_SCS_GENERALIZED_DIAGONAL,
    },
    tuple(['-wd', '--workdir']): {
        'help': phelp_scs_workdir,
        'type': str,
        'default': DEFAULT_SCS_GENERALIZED_WORKDIR,
    },
    tuple(['--seed']): {
        'help': phelp_scs_seed,
        'type': int,
        'default': None,
    },
}

SCSGeneralized_action_args_dict = {**action_args_dict}
