from .common import *
#
from .SignedGraph import *
#
from .phelp.SignedGraph import *
from .phelp.Lattice3D import *
from .phelp.TransCluster import *
#
from .defs.SignedGraph import *
from .defs.Lattice3D import *
from .defs.TransCluster import *
#
L3D_args = {
    'L': {
        'help': phelp_L,
        'type': int,
        'nargs': '+',
        'action': IntOrTriple
    },
    **SG_args
}
L3D_opt_args = {
    tuple(['-g', '--geometry']): {
        'help': phelp_geo,
        'type': str,
        'default': DEFAULT_GEO
    },
    tuple(['--pdil']): {
        'help': phelp_pdil,
        'type': float,
        'default': DEFAULT_L3D_TRANSCLUSTER_PDIL,
    },
    tuple(['--mu']): {
        'help': phelp_mu,
        'type': float,
        'default': DEFAULT_L3D_TRANSCLUSTER_MU,
    },
    tuple(['--sigma']): {
        'help': phelp_sigma,
        'type': float,
        'default': DEFAULT_L3D_TRANSCLUSTER_SIGMA,
    },
    tuple(['--edge_weight']): {
        'help': phelp_edge_weight,
        'type': str,
        'default': DEFAULT_L3D_TRANSCLUSTER_EDGE_WEIGHT,
    },
    tuple(['-na', '--number_of_averages']): {
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_L3D_TRANSCLUSTER_NAVG,
    },
    **SG_opt_args,
}