from .common import *
from .phelp.TransCluster import *
from .defs.TransCluster import *
from .Lattice2D import *

L2D_TransCluster_progName = 'L2D_TransCluster'
L2D_TransCluster_progNameShrt = 'L2DTC'
L2D_TransCluster_description = """
    Transient cluster analysis utilities for 2D signed lattices: {L2D_TransCluster_progName}.py
"""

L2D_TransCluster_args = {**L2D_args}

L2D_TransCluster_optional_args_dict = {
    tuple(['--prew']): {
        'help': phelp_prew,
        'type': float,
        'default': DEFAULT_L2D_TRANSCLUSTER_PREW,
    },
    tuple(['-m', '--mode']): {
        'help': phelp_transcluster_mode,
        'type': str,
        'default': DEFAULT_L2D_TRANSCLUSTER_MODE,
    },
    tuple(['-sf', '--save_frequency']): {
        'help': phelp_data_save_freq,
        'type': int,
        'default': DEFAULT_L2D_TRANSCLUSTER_SAVE_FREQUENCY,
    },
    tuple(['-o', '--out_suffix']): {
        'help': phelp_out_suffix,
        'type': str,
        'default': DEFAULT_L2D_TRANSCLUSTER_OUT_SUFFIX,
    },
    tuple(['-t', '--float_type']): {
        'help': phelp_float_type,
        'type': str,
        'default': DEFAULT_L2D_TRANSCLUSTER_FLOAT_TYPE,
    },
}

L2D_TransCluster_action_args_dict = {**action_args_dict}
