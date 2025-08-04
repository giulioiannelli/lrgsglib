from .common import *

from .phelp.SignedGraph import *
from .defs.SignedGraph import *

SG_args = {
    'p': {
        'help': phelp_p,
        'type': float
    }
}

SG_opt_args = {
    tuple(['-c', '--cell_type']): {
        'help': phelp_cell,
        'type': str,
        'default': DEFAULT_CELL
    },
    tuple(['-na', '--number_of_averages']): {
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_NAVG
    },
    tuple(['-wd', '--workdir']): {
        'help': phelp_workdir,
        'type': str,
        'default': DEFAULT_WORKDIR
    },
}