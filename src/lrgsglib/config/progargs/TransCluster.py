from .common import *
from .phelp.TransCluster import *
from .defs.TransCluster import *
from .Lattice2D import *
from .SCSGeneralized import *

L2D_TransCluster_progName = 'L2D_TransCluster'
L2D_TransCluster_progNameShrt = 'L2DTC'
L2D_TransCluster_description = f"""
    Transition and cluster analysis utilities for 2D signed lattices: {L2D_TransCluster_progName}.py
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

SCS_TransCluster_progName = 'SCS_TransCluster'
SCS_TransCluster_progNameShrt = 'SCSTC'
SCS_TransCluster_description = f"""
    Order parameter and clustering analysis for generalized SCS networks: {SCS_TransCluster_progName}.py
"""

SCS_TransCluster_args = {**SCSGeneralized_args}

SCS_TransCluster_optional_args_dict = {
    tuple(['--J0']): {
        'help': phelp_scs_j0,
        'type': float,
        'default': DEFAULT_SCS_GENERALIZED_J0,
    },
    tuple(['-na', '--number_of_averages']): {
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_SCS_TRANSCLUSTER_NAVG,
    },
    tuple(['-m', '--mode']): {
        'help': phelp_transcluster_mode,
        'type': str,
        'default': DEFAULT_SCS_TRANSCLUSTER_MODE,
    },
    tuple(['-sf', '--save_frequency']): {
        'help': phelp_scs_save_frequency,
        'type': int,
        'default': DEFAULT_SCS_TRANSCLUSTER_SAVE_FREQUENCY,
    },
    tuple(['-o', '--out_suffix']): {
        'help': phelp_out_suffix,
        'type': str,
        'default': DEFAULT_SCS_TRANSCLUSTER_OUT_SUFFIX,
    },
    tuple(['-t', '--float_type']): {
        'help': phelp_scs_float_type,
        'type': str,
        'default': DEFAULT_SCS_TRANSCLUSTER_FLOAT_TYPE,
    },
    tuple(['--backend']): {
        'help': phelp_scs_backend,
        'type': str,
        'default': DEFAULT_SCS_TRANSCLUSTER_BACKEND,
    },
    tuple(['--partition-rule']): {
        'help': phelp_scs_partition_rule,
        'type': str,
        'default': DEFAULT_SCS_TRANSCLUSTER_PARTITION_RULE,
    },
}

SCS_TransCluster_action_args_dict = {**SCSGeneralized_action_args_dict}
