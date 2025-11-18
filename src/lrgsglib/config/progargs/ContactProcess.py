import argparse

from .common import *
from .Lattice2D import *

from .phelp.ContactProcess import *
from .defs.ContactProcess import *

Contact_args = {
    tuple(['-dy', '--dynamics']): {
        'help': phelp_dynamics,
        'type': str,
        'choices': ["EI", "SIR"],
        'default': DEFAULT_DYNAMICS,
    },
}

Contact_opt_args = {
    tuple(['-mu', '--mu']): {
        'help': phelp_mu,
        'type': float,
        'default': DEFAULT_MU,
    },
    tuple(['-ga', '--gamma']): {
        'help': phelp_gamma,
        'type': float,
        'default': DEFAULT_GAMMA,
    },
    tuple(['-ac', '--activation']): {
        'help': phelp_activation,
        'type': str,
        'default': DEFAULT_ACTIVATION,
    },
    tuple(['--save_density']): {
        'help': phelp_save_density,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SAVE_DENSITY,
    },
    tuple(['-st', '--state_type']): {
        'help': phelp_state_type,
        'type': str,
        'default': DEFAULT_STATE_TYPE,
    },
    tuple(['-ns', '--num_log_samples']): {
        'help': phelp_num_log_samples,
        'type': int,
        'default': DEFAULT_NUM_LOG_SAMPLES,
    },
    tuple(['-sp', '--simpref']): {
        'help': phelp_simpref,
        'type': int,
        'default': DEFAULT_SIM_PREF,
    },
    tuple(['-ic', '--init_cond']): {
        'help': phelp_ic_contact,
        'type': str,
        'default': DEFAULT_INIT_COND,
    },
    tuple(['-rl', '--runlang']): {
        'help': phelp_runlang_contact,
        'type': str,
        'default': DEFAULT_RUNLANG,
    },
    tuple(['-os', '--out_suffix']): {
        'help': phelp_outsuffix_contact,
        'type': str,
        'default': DEFAULT_OUTSFFX,
    },
    tuple(['-rnds', '--randstr']): {
        'help': phelp_randstr_contact,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_randstr,
    },
}

L2D_CPROC_progname = "L2D_ContactProcess"
L2D_CPROC_progname_shrt = "L2DCP"
L2D_CPROC_description = f"""
    Computational resources for running the ContactProcess dynamics on
    2D lattices: {L2D_CPROC_progname}.py
"""

L2D_CPROC_args = {**L2D_args, **Contact_args}
L2D_CPROC_opt_args = {**L2D_opt_args, **Contact_opt_args}
L2D_CPROC_action_args = {**action_args_dict}
