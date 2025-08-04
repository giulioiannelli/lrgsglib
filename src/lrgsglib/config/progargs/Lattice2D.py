from .common import *
#
from .phelp.SignedGraph import *
from .phelp.Lattice2D import *
#
from .defs.SignedGraph import *
from .defs.Lattice2D import *
#
L2D_args = {
    'L': {
        'help': phelp_L, 
        'type': int
    },
    'p': {
        'help': phelp_p,
        'type': float
    }
}
L2D_opt_args = {
    tuple(['-c', '--cell_type']): {
        'help': phelp_cell,
        'type': str,
        'default': DEFAULT_CELL
    },
    tuple(['-g', '--geometry']): {
        'help': phelp_geo,
        'type': str,
        'default': DEFAULT_GEO
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
    tuple(['-cpt', '--compute']): {
        'help': phelp_compute,
        'type': str,
        'default': DEFAULT_COMPUTE
    },
}