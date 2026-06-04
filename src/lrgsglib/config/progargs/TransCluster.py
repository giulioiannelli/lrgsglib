from .phelp.TransCluster import *
from .defs.TransCluster import *
from .Lattice2D import *
from .Lattice3D import *
from .SCSNN import *

# Laplacian type selector (canonical source: config/const.py; help string reused
# from the SlaplSpect program so the wording stays in one place).
from ..const import SG_LAPL_DEFAULT_TYPE, SG_LAPL_TYPES
from .phelp.SlaplSpect import phelp_laplacian_type

def _transcluster_common_optional_args(
    *,
    mode_default,
    save_frequency_default,
    out_suffix_default,
    float_type_default,
    save_frequency_help,
    out_suffix_help=phelp_out_suffix,
    float_type_help=phelp_float_type,
):
    return {
        tuple(["-m", "--mode"]): {
            "help": phelp_transcluster_mode,
            "type": str,
            "default": mode_default,
        },
        tuple(["-sf", "--save_frequency"]): {
            "help": save_frequency_help,
            "type": int,
            "default": save_frequency_default,
        },
        tuple(["-o", "--out_suffix"]): {
            "help": out_suffix_help,
            "type": str,
            "default": out_suffix_default,
        },
        tuple(["-t", "--float_type"]): {
            "help": float_type_help,
            "type": str,
            "default": float_type_default,
        },
    }


L2D_TransCluster_progName = 'L2D_TransCluster'
L2D_TransCluster_progNameShrt = 'L2DTC'
L2D_TransCluster_description = f"""
    Phase Transition and cluster analysis utilities for 2D signed lattices:
    {L2D_TransCluster_progName}.py
"""
L2D_TransCluster_args = {**L2D_args}
L2D_TransCluster_optional_args_dict = {
    **L2D_opt_args,
    **_transcluster_common_optional_args(
        mode_default=DEFAULT_L2D_TRANSCLUSTER_MODE,
        save_frequency_default=DEFAULT_L2D_TRANSCLUSTER_SAVE_FREQUENCY,
        out_suffix_default=DEFAULT_L2D_TRANSCLUSTER_OUT_SUFFIX,
        float_type_default=DEFAULT_L2D_TRANSCLUSTER_FLOAT_TYPE,
        save_frequency_help=phelp_data_save_freq,
    ),
    tuple(["-lt", "--laplacian-type"]): {
        "help": phelp_laplacian_type,
        "type": str,
        "default": SG_LAPL_DEFAULT_TYPE,
        "choices": list(SG_LAPL_TYPES),
    },
}
L2D_TransCluster_action_args_dict = {**action_args_dict}

L2D_TransCluster_srun_description = (
    f"Serializer for {L2D_TransCluster_progName}.py"
)

L2D_TransCluster_srun_optional_args_dict = {
    tuple(["-m", "--mode"]): {
        "help": (
            "Execution mode: 'pCluster', 'ordParam', 'slanzarv_pCluster', "
            "or 'slanzarv_ordParam'. If 'slanzarv_' prefix is present, "
            "jobs will be submitted via slanzarv."
        ),
        "type": str,
        "choices": [
            "pCluster",
            "ordParam",
            "slanzarv_pCluster",
            "slanzarv_ordParam",
        ],
        "default": "slanzarv_ordParam",
    },
    tuple(["--L-list"]): {
        "help": (
            "List of lattice sizes | "
            f"default={list(DEFAULT_L2D_TRANSCLUSTER_SRUN_L_LIST)}"
        ),
        "type": int,
        "nargs": "+",
        "default": None,
    },
    tuple(["--L-linsp"]): {
        "help": "Semicolon-separated linspace tuples for lattice sizes",
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--p-list"]): {
        "help": "Explicit probabilities for flipped edges",
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--p-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for flipped-edge "
            "probabilities"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--mu-list"]): {
        "help": (
            "Explicit mu values for normal edge weights "
            "(triggers normal mode)"
        ),
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--mu-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for mu values "
            "(triggers normal mode)"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--sigma-list"]): {
        "help": (
            "Explicit sigma values for normal edge weights "
            "(triggers normal mode)"
        ),
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--sigma-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for sigma values "
            "(triggers normal mode)"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    **srun_opt_args,
}
#
L3D_TransCluster_progName = 'L3D_TransCluster'
L3D_TransCluster_progNameShrt = 'L3DTC'
L3D_TransCluster_description = f"""
    Phase Transition and cluster analysis utilities for 3D signed lattices:
    {L3D_TransCluster_progName}.py
"""

L3D_TransCluster_args = {**L3D_args}

L3D_TransCluster_optional_args_dict = {
    **L3D_opt_args,
    **_transcluster_common_optional_args(
        mode_default=DEFAULT_L3D_TRANSCLUSTER_MODE,
        save_frequency_default=DEFAULT_L3D_TRANSCLUSTER_SAVE_FREQUENCY,
        out_suffix_default=DEFAULT_L3D_TRANSCLUSTER_OUT_SUFFIX,
        float_type_default=DEFAULT_L3D_TRANSCLUSTER_FLOAT_TYPE,
        save_frequency_help=phelp_data_save_freq,
    ),
}

L3D_TransCluster_action_args_dict = {**action_args_dict}

L3D_TransCluster_srun_description = (
    f"Serializer for {L3D_TransCluster_progName}.py"
)

L3D_TransCluster_srun_optional_args_dict = {
    tuple(["-m", "--mode"]): {
        "help": (
            "Execution mode: 'pCluster', 'ordParam', 'slanzarv_pCluster', "
            "or 'slanzarv_ordParam'. If 'slanzarv_' prefix is present, "
            "jobs will be submitted via slanzarv."
        ),
        "type": str,
        "choices": [
            "pCluster",
            "ordParam",
            "slanzarv_pCluster",
            "slanzarv_ordParam",
        ],
        "default": "slanzarv_ordParam",
    },
    tuple(["--L-list"]): {
        "help": (
            "List of lattice sizes | "
            f"default={list(DEFAULT_L3D_TRANSCLUSTER_SRUN_L_LIST)}"
        ),
        "type": int,
        "nargs": "+",
        "default": None,
    },
    tuple(["--L-linsp"]): {
        "help": "Semicolon-separated linspace tuples for lattice sizes",
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--p-list"]): {
        "help": "Explicit probabilities for flipped edges",
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--p-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for flipped-edge "
            "probabilities"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--mu-list"]): {
        "help": (
            "Explicit mu values for normal edge weights "
            "(triggers normal mode)"
        ),
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--mu-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for mu values "
            "(triggers normal mode)"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--sigma-list"]): {
        "help": (
            "Explicit sigma values for normal edge weights "
            "(triggers normal mode)"
        ),
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--sigma-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for sigma values "
            "(triggers normal mode)"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    **srun_opt_args,
}

SCS_TransCluster_progName = 'SCS_TransCluster'
SCS_TransCluster_progNameShrt = 'SCSTC'
SCS_TransCluster_description = f"""
    Order parameter and clustering analysis for generalized SCS networks:
    {SCS_TransCluster_progName}.py
"""

SCS_TransCluster_args = {**SCSGeneralized_args}

SCS_TransCluster_optional_args_dict = {
    tuple(['--J0']): {
        'help': phelp_scs_j0,
        'type': float,
        'default': DEFAULT_SCS_NN_J0,
    },
    tuple(['-na', '--number_of_averages']): {
        'help': phelp_navg,
        'type': int,
        'default': DEFAULT_SCS_NN_NAVG,
    },
    **_transcluster_common_optional_args(
        mode_default=DEFAULT_SCS_NN_MODE,
        save_frequency_default=DEFAULT_SCS_NN_SAVE_FREQUENCY,
        out_suffix_default=DEFAULT_SCS_NN_OUT_SUFFIX,
        float_type_default=DEFAULT_SCS_NN_FLOAT_TYPE,
        save_frequency_help=phelp_scs_save_frequency,
        float_type_help=phelp_scs_float_type,
    ),
    tuple(['--backend']): {
        'help': phelp_scs_backend,
        'type': str,
        'default': DEFAULT_SCS_NN_BACKEND,
    },
    tuple(['--partition-rule']): {
        'help': phelp_scs_partition_rule,
        'type': str,
        'default': DEFAULT_SCS_NN_PARTITION_RULE,
    },
}

SCS_TransCluster_action_args_dict = {**SCSGeneralized_action_args_dict}

SCS_TransCluster_srun_description = (
    f"Serializer for {SCS_TransCluster_progName}.py"
)

SCS_TransCluster_srun_optional_args_dict = {
    tuple(["-m", "--mode"]): {
        "help": (
            "Execution mode: 'pCluster', 'ordParam', 'slanzarv_pCluster', "
            "or 'slanzarv_ordParam'. If 'slanzarv_' prefix is present, "
            "jobs will be submitted via slanzarv."
        ),
        "type": str,
        "choices": [
            "pCluster",
            "ordParam",
            "slanzarv_pCluster",
            "slanzarv_ordParam",
        ],
        "default": "slanzarv_ordParam",
    },
    tuple(["--N-list"]): {
        "help": "List of network sizes N",
        "type": int,
        "nargs": "+",
        "default": list(DEFAULT_SCS_NN_SRUN_N_LIST),
    },
    tuple(["--gamma-list"]): {
        "help": "List of gamma values",
        "type": float,
        "nargs": "+",
        "default": list(DEFAULT_SCS_NN_SRUN_GAMMA_LIST),
    },
    tuple(["--J0-list"]): {
        "help": "List of J0 values",
        "type": float,
        "nargs": "+",
        "default": list(DEFAULT_SCS_NN_SRUN_J0_LIST),
    },
    tuple(["--J0-linsp"]): {
        "help": "Semicolon-separated linspace tuples for J0 values",
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--J-list"]): {
        "help": "List of J coupling values",
        "type": float,
        "nargs": "+",
        "default": list(DEFAULT_SCS_NN_SRUN_J_LIST),
    },
    tuple(["--J-linsp"]): {
        "help": "Semicolon-separated linspace tuples for J values",
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--g-list"]): {
        "help": "List of g coupling values",
        "type": float,
        "nargs": "+",
        "default": list(DEFAULT_SCS_NN_SRUN_G_LIST),
    },
    tuple(["--g-linsp"]): {
        "help": "Semicolon-separated linspace tuples for g values",
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--diagonal-list"]): {
        "help": "List of diagonal values",
        "type": str,
        "nargs": "+",
        "default": list(DEFAULT_SCS_NN_SRUN_DIAGONAL_LIST),
    },
    **srun_opt_args,
}
