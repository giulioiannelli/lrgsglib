"""argparse dicts for the SignedRW walker CLI (``L2D_SignedRW.py``).

Exposes ``SRW_opt_args`` (walker-specific options) and the composite
``L2D_SRW_*`` bundles used by ``parsers/L2D_SignedRW.py`` and the
downstream serializer.
"""

import argparse

from .common import *
from .Lattice2D import *
from .phelp.srw import *
from .defs.srw import *


SRW_opt_args = {
    tuple(['-ru', '--rule']): {
        'help': phelp_srw_rule,
        'type': str,
        'choices': list(SRW_RULES),
        'default': DEFAULT_SRW_RULE,
    },
    tuple(['-m', '--mode']): {
        'help': phelp_srw_mode,
        'type': str,
        'choices': list(SRW_MODES),
        'default': DEFAULT_SRW_MODE,
    },
    tuple(['-nw', '--n_walkers']): {
        'help': phelp_srw_n_walkers,
        'type': int,
        'default': DEFAULT_SRW_N_WALKERS,
    },
    tuple(['-cv', '--coverage']): {
        'help': phelp_srw_coverage,
        'type': float,
        'default': DEFAULT_SRW_COVERAGE_FRAC,
    },
    tuple(['-xn', '--x_node']): {
        'help': phelp_srw_x_node,
        'type': str,
        'choices': list(SRW_X_NODE_BEHAVIORS),
        'default': DEFAULT_SRW_X_NODE,
    },
    tuple(['--max_n_cross']): {
        'help': phelp_srw_max_n_cross,
        'type': int,
        'default': DEFAULT_SRW_MAX_N_CROSS,
    },
    tuple(['--start_a']): {
        'help': phelp_srw_start_a,
        'type': str,
        'choices': list(SRW_START_PROTOCOLS),
        'default': DEFAULT_SRW_START,
    },
    tuple(['--start_b']): {
        'help': phelp_srw_start_b,
        'type': str,
        'choices': list(SRW_START_PROTOCOLS),
        'default': DEFAULT_SRW_START_B,
    },
    tuple(['--start_a_node']): {
        'help': phelp_srw_start_node_a,
        'type': int,
        'default': None,
    },
    tuple(['--start_b_node']): {
        'help': phelp_srw_start_node_b,
        'type': int,
        'default': None,
    },
    tuple(['-sa', '--seed_a']): {
        'help': phelp_srw_seed_a,
        'type': int,
        'default': DEFAULT_SRW_SEED,
    },
    tuple(['-sb', '--seed_b']): {
        'help': phelp_srw_seed_b,
        'type': int,
        'default': DEFAULT_SRW_SEED_B,
    },
    tuple(['--store_trajectory']): {
        'help': phelp_srw_store_trajectory,
        'action': argparse.BooleanOptionalAction,
        'default': False,
    },
    tuple(['--store_per_walker_visits']): {
        'help': phelp_srw_store_per_walker_visits,
        'action': argparse.BooleanOptionalAction,
        'default': False,
    },
    tuple(['-rl', '--runlang']): {
        'help': phelp_srw_runlang,
        'type': str,
        'default': 'py',
    },
}


L2D_SRW_progname = "L2D_SignedRW"
L2D_SRW_progname_shrt = "L2DSRW"
L2D_SRW_description = f"""
    Signed random-walker dynamics on 2D lattices: {L2D_SRW_progname}.py.
    Runs an ensemble of single-walker simulations (absorb / kill / sticky
    rule) and, in --mode=pair, computes the visit-count overlap between
    two independent ensembles with distinct start configurations.
"""

L2D_SRW_args = {**L2D_args}
L2D_SRW_opt_args = {**L2D_opt_args, **SRW_opt_args}
L2D_SRW_action_args = {**action_args_dict}
