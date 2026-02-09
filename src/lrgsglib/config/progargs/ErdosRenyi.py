from .common import *
#
from .phelp.SignedGraph import *
from .phelp.ErdosRenyi import *
#
from .defs.SignedGraph import *
from .defs.ErdosRenyi import *
#
ER_args = {
    'N': {
        'help': phelp_N,
        'type': int
    },
    'p': {
        'help': phelp_p_er,
        'type': float
    },
    'pflip': {
        'help': phelp_pflip,
        'type': float
    }
}
ER_opt_args = {
    tuple(['-c', '--cell_type']): {
        'help': phelp_cell,
        'type': str,
        'default': DEFAULT_CELL
    },
    tuple(['-na', '--number_of_averages']): {
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_NAVG_ER
    },
    tuple(['-wd', '--workdir']): {
        'help': phelp_workdir,
        'type': str,
        'default': DEFAULT_WORKDIR
    },
    tuple(['-ge', '--graph_engine']): {
        'help': phelp_graph_engine,
        'type': str,
        'default': DEFAULT_GRAPH_ENGINE,
        'choices': ['nx', 'gt'],
    },
}
