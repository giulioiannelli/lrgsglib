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
    tuple(['-sf', '--save_frequency']): {
        'help': phelp_save_frequency,
        'type': int,
        'default': DEFAULT_SAVE_FREQUENCY,
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
    tuple(['-esdt', '--early_stop_density_threshold']): {
        'help': 'Early stopping density threshold. If the average of the last M timesteps (M=10%% of time series length) exceeds this value after 20 runs, stop the simulation.',
        'type': float,
        'default': None,
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

L2D_ContactProcess_srun_description = (
    f"Serialiser for {L2D_CPROC_progname}.py"
)

L2D_ContactProcess_srun_optional_args_dict = {
    tuple(["-m", "--mode"]): {
        "help": (
            "Execution mode; prefix with 'slanzarv_' to submit via slanzarv "
            "(e.g. slanzarv_EI)."
        ),
        "type": str,
        "choices": ["EI", "SIR", "slanzarv_EI", "slanzarv_SIR"],
        "default": "slanzarv_EI",
    },
    tuple(["--L-list"]): {
        "help": (
            "List of lattice sizes | "
            f"default={list(DEFAULT_L2D_CONTACTPROCESS_SRUN_L_LIST)}"
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
        "help": (
            "Contact process occupation probabilities | "
            f"default={list(DEFAULT_L2D_CONTACTPROCESS_SRUN_P_LIST)}"
        ),
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--p-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for occupation probabilities"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--gamma-list"]): {
        "help": (
            "Coupling values for EI dynamics | "
            f"default={list(DEFAULT_L2D_CONTACTPROCESS_SRUN_GAMMA_LIST)}"
        ),
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--gamma-linsp"]): {
        "help": "Semicolon-separated linspace tuples for EI coupling values",
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["--mu-list"]): {
        "help": (
            "Infection rates μ for SIR dynamics | "
            f"default={list(DEFAULT_L2D_CONTACTPROCESS_SRUN_MU_LIST)}"
        ),
        "type": float,
        "nargs": "+",
        "default": None,
    },
    tuple(["--mu-linsp"]): {
        "help": (
            "Semicolon-separated linspace tuples for infection-rate μ values"
        ),
        "type": parse_multiple_linspace,
        "default": None,
    },
    tuple(["-rl", "--runlang"]): {
        "help": phelp_runlang_contact,
        "type": str,
        "default": DEFAULT_RUNLANG,
    },
    **srun_opt_args,
}
