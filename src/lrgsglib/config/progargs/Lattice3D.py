from .common import *
#
from .SignedGraph import *
#
from .phelp.SignedGraph import *
from .phelp.Lattice3D import *
#
from .defs.SignedGraph import *
from .defs.Lattice3D import *
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
    **SG_opt_args,
}