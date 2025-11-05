from __future__ import annotations

import argparse
from typing import Sequence

from lrgsglib.config.progargs import *

DEFAULT_SCS_TC_SERIALISER_N_LIST: Sequence[int] = (500,)
DEFAULT_SCS_TC_SERIALISER_GAMMA_LIST: Sequence[float] = (1.0,)
DEFAULT_SCS_TC_SERIALISER_J0_LIST: Sequence[float] = (DEFAULT_SCS_NN_J0,)

_description = f"""
    Serialiser for {SCS_TransCluster_progName}.py
"""

_list_args = {
    ("--N-list",): {
        "help": f"List of network sizes N | default={list(DEFAULT_SCS_TC_SERIALISER_N_LIST)}",
        "type": int,
        "nargs": "+",
        "default": DEFAULT_SCS_TC_SERIALISER_N_LIST,
    },
    ("--gamma-list",): {
        "help": f"List of gamma values | default={list(DEFAULT_SCS_TC_SERIALISER_GAMMA_LIST)}",
        "type": float,
        "nargs": "+",
        "default": DEFAULT_SCS_TC_SERIALISER_GAMMA_LIST,
    },
    ("--J0-list",): {
        "help": f"List of J0 values | default={list(DEFAULT_SCS_TC_SERIALISER_J0_LIST)}",
        "type": float,
        "nargs": "+",
        "default": DEFAULT_SCS_TC_SERIALISER_J0_LIST,
    },
}

_scalar_args = {
    ("--J",): {
        "help": phelp_scs_J,
        "type": float,
        "default": DEFAULT_SCS_NN_J,
    },
    ("--g",): {
        "help": phelp_scs_g,
        "type": float,
        "default": DEFAULT_SCS_NN_G,
    },
    ("--diagonal",): {
        "help": phelp_scs_diagonal,
        "type": str,
        "default": DEFAULT_SCS_NN_DIAGONAL,
    },
    ("-wd", "--workdir"): {
        "help": phelp_scs_workdir,
        "type": str,
        "default": DEFAULT_SCS_NN_WORKDIR,
    },
    ("-na", "--number_of_averages"): {
        "help": phelp_navg,
        "type": int,
        "default": DEFAULT_SCS_TRANSCLUSTER_NAVG,
    },
    ("-o", "--out_suffix"): {
        "help": phelp_out_suffix,
        "type": str,
        "default": DEFAULT_SCS_TRANSCLUSTER_OUT_SUFFIX,
    },
    ("-m", "--mode"): {
        "help": phelp_transcluster_mode,
        "type": str,
        "default": DEFAULT_SCS_TRANSCLUSTER_MODE,
    },
    ("-sf", "--save_frequency"): {
        "help": phelp_scs_save_frequency,
        "type": int,
        "default": DEFAULT_SCS_TRANSCLUSTER_SAVE_FREQUENCY,
    },
    ("-t", "--float_type"): {
        "help": phelp_scs_float_type,
        "type": str,
        "default": DEFAULT_SCS_TRANSCLUSTER_FLOAT_TYPE,
    },
    ("--backend",): {
        "help": phelp_scs_backend,
        "type": str,
        "default": DEFAULT_SCS_TRANSCLUSTER_BACKEND,
    },
    ("--partition-rule",): {
        "help": phelp_scs_partition_rule,
        "type": str,
        "default": DEFAULT_SCS_TRANSCLUSTER_PARTITION_RULE,
    },
    ("--seed",): {
        "help": phelp_scs_seed,
        "type": int,
        "default": None,
    },
}

_parser = argparse.ArgumentParser(
    description=_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
)

for key, opts in {**_list_args, **_scalar_args}.items():
    _parser.add_argument(*key, **opts)

for key, opts in {**srun_opt_args, **srun_action_args}.items():
    _parser.add_argument(*key, **opts)

parser = _parser
