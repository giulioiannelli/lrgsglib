from .common import *

from .phelp.SignedGraph import *
from .defs.SignedGraph import *

SG_args: dict = {}

SG_opt_args = {
    tuple(['-p', '--pflip']): {
        'help': phelp_p,
        'type': float,
        'default': 0.0,
        'dest': 'p',
    },
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