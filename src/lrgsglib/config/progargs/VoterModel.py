import argparse

from .common import *
from .Lattice2D import *

# Shared serializer infrastructure (sweep / dispatch / slanzarv / save help and
# the generic array-dispatch defaults) is reused from the IsingDynamics program
# args -- it is not Ising-specific. Reuse, do not duplicate (hard-ban #1).
from .phelp.IsingDynamics import *
from .defs.IsingDynamics import (
    DISPATCH_CHOICES,
    DEFAULT_DISPATCH,
    DEFAULT_ARRAY_CONCURRENT,
    DEFAULT_ARRAY_PARTITION,
    DEFAULT_ARRAY_TIME,
    DEFAULT_QUENCH_ID,
)

from .phelp.VoterModel import *
from .defs.VoterModel import *

# Voter dynamics defaults / enumerations are single-sourced in the model
# defaults module (pure constants, no lrgsglib imports -> no cycle).
from lrgsglib.statsys.VoterModel.defaults import (
    VOTER_RULES,
    DEFAULT_RULE,
    VOTER_UPD_MODES,
    DEFAULT_UPD_MODE,
    DEFAULT_QVOTER_Q,
    DEFAULT_NOISE_EPS,
    DEFAULT_NONLIN_ALPHA,
    CLUSTER_MODES,
    DEFAULT_CLUSTER_MODE,
    DEFAULT_ABSORBING_EVERY,
    DEFAULT_SOUT_EVERY,
)

# ===========================================================================
# Program args: L2D_VoterModel.py
# ===========================================================================
# No required positional beyond the lattice side ``L`` (from L2D_args).
Voter_args: dict = {}

# Voter dynamics + run options (geometry comes from L2D_opt_args at the parser
# level, exactly like IsingDynamics).
Voter_opt_args = {
    # --- Axis A: local update rule ---
    tuple(['--rule']): {
        'help': phelp_voter_rule,
        'type': str,
        'choices': list(VOTER_RULES),
        'default': DEFAULT_RULE,
    },
    tuple(['-q', '--q']): {
        'help': phelp_voter_q,
        'type': int,
        'default': DEFAULT_QVOTER_Q,
    },
    tuple(['--eps']): {
        'help': phelp_voter_eps,
        'type': float,
        'default': DEFAULT_NOISE_EPS,
    },
    tuple(['--alpha']): {
        'help': phelp_voter_alpha,
        'type': float,
        'default': DEFAULT_NONLIN_ALPHA,
    },
    # --- Axis B: update schedule ---
    tuple(['-upd', '--upd_mode']): {
        'help': phelp_voter_upd_mode,
        'type': str,
        'choices': list(VOTER_UPD_MODES),
        'default': DEFAULT_UPD_MODE,
    },
    # --- backend / init / time ---
    tuple(['-rl', '--runlang']): {
        'help': phelp_voter_runlang,
        'type': str,
        'default': DEFAULT_VOTER_RUNLANG,
    },
    tuple(['-ic', '--init_cond']): {
        'help': phelp_voter_ic,
        'type': str,
        'default': DEFAULT_VOTER_INIT_COND,
    },
    tuple(['--steps', '-stp']): {
        'help': phelp_voter_steps,
        'type': int,
        'default': DEFAULT_VOTER_STEPS,
    },
    tuple(['--simref', '-sp']): {
        'help': phelp_voter_simref,
        'type': float,
        'default': DEFAULT_VOTER_SIMREF,
    },
    tuple(['-st', '--state_type']): {
        'help': phelp_voter_state_type,
        'type': str,
        'default': DEFAULT_VOTER_STATE_TYPE,
    },
    tuple(['-fq', '--freq']): {
        'help': phelp_voter_freq,
        'type': int,
        'default': DEFAULT_VOTER_FREQ,
    },
    tuple(['-ns', '--num_log_samples']): {
        'help': phelp_voter_num_log_samples,
        'type': int,
        'default': DEFAULT_VOTER_NUM_LOG_SAMPLES,
    },
    tuple(['--seed']): {
        'help': phelp_voter_seed,
        'type': int,
        'default': DEFAULT_VOTER_SEED,
    },
    # --- absorbing-state detection ---
    tuple(['--absorbing-check']): {
        'help': phelp_voter_absorbing_check,
        'action': argparse.BooleanOptionalAction,
        'default': None,
    },
    tuple(['--absorbing-every']): {
        'help': phelp_voter_absorbing_every,
        'type': int,
        'default': DEFAULT_ABSORBING_EVERY,
    },
    # --- cluster-size distribution ---
    tuple(['-cm', '--cluster-mode']): {
        'help': phelp_voter_cluster_mode,
        'type': str,
        'choices': list(CLUSTER_MODES),
        'default': DEFAULT_CLUSTER_MODE,
    },
    # --- trajectory subsampling ---
    tuple(['--sout-every']): {
        'help': phelp_voter_sout_every,
        'type': int,
        'default': DEFAULT_SOUT_EVERY,
    },
    tuple(['--sout-nlog']): {
        'help': phelp_voter_sout_nlog,
        'type': int,
        'default': DEFAULT_VOTER_SOUT_NLOG,
    },
    tuple(['--sout-force-full']): {
        'help': phelp_voter_sout_force_full,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VOTER_SOUT_FORCE_FULL,
    },
    # --- single-quench (SLURM array) mode ---
    tuple(['-qid', '--quench-id']): {
        'help': phelp_quench_id,
        'type': int,
        'default': DEFAULT_QUENCH_ID,
    },
    # --- output / persistence ---
    tuple(['-os', '--out_suffix']): {
        'help': phelp_outsuffix,
        'type': str,
        'default': DEFAULT_OUTSFFX,
    },
    tuple(['-rnds', '--randstr']): {
        'help': phelp_voter_randstr,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VOTER_RANDSTR,
    },
}

# Observable save gates (one per VoterModel ObservableSet member).
Voter_save_args = {
    tuple(['--save-magn']): {
        'help': phelp_voter_save_magn,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VOTER_SAVE_MAGN,
    },
    tuple(['--save-sout']): {
        'help': phelp_voter_save_sout,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VOTER_SAVE_SOUT,
    },
    tuple(['--save-cldist']): {
        'help': phelp_voter_save_cldist,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VOTER_SAVE_CLDIST,
    },
    tuple(['--save-all']): {
        'help': phelp_voter_save_all,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VOTER_SAVE_ALL,
    },
}

# ===========================================================================
# Serializer args: L2D_VoterModel_Serializer.py
# ===========================================================================
# Sweep axes. ``side1_list`` and ``pflip_linsp`` are always swept (reusing the
# shared lattice defaults); the voter dynamics axes are opt-in (default None =
# forward the single scalar from the program flags, no extra nesting).
Voter_srun_list_args = {
    tuple(['-s1', '--side1_list']): {
        'help': phelp_side1_list,
        'type': int,
        'nargs': '+',
        'default': DEFAULT_SIDE1_LIST,
    },
    tuple(['-pFT', '--pflip_linsp']): {
        'help': phelp_pflip_linsp,
        'type': parse_multiple_linspace,
        'default': DEFAULT_PFLIP_LINSP,
    },
    tuple(['--rule-list']): {
        'help': phelp_voter_rule_list,
        'type': str,
        'nargs': '+',
        'choices': list(VOTER_RULES),
        'default': None,
    },
    tuple(['--upd-mode-list']): {
        'help': phelp_voter_upd_mode_list,
        'type': str,
        'nargs': '+',
        'choices': list(VOTER_UPD_MODES),
        'default': None,
    },
    tuple(['--alpha-list']): {
        'help': phelp_voter_alpha_list,
        'type': float,
        'nargs': '+',
        'default': None,
    },
    tuple(['--alpha-linsp']): {
        'help': phelp_voter_alpha_linsp,
        'type': parse_multiple_linspace,
        'default': None,
    },
    tuple(['--eps-list']): {
        'help': phelp_voter_eps_list,
        'type': float,
        'nargs': '+',
        'default': None,
    },
    tuple(['--eps-linsp']): {
        'help': phelp_voter_eps_linsp,
        'type': parse_multiple_linspace,
        'default': None,
    },
    tuple(['--q-list']): {
        'help': phelp_voter_q_list,
        'type': int,
        'nargs': '+',
        'default': None,
    },
}

# Array-dispatch (SLURM --array) controls; mirror IsingDynamics but without the
# Ising-specific -mL/-cL list axes.
Voter_srun_array_args = {
    tuple(['--dispatch']): {
        'help': phelp_dispatch,
        'type': str,
        'choices': DISPATCH_CHOICES,
        'default': DEFAULT_DISPATCH,
    },
    tuple(['-ac', '--array-concurrent']): {
        'help': phelp_array_concurrent,
        'type': int,
        'default': DEFAULT_ARRAY_CONCURRENT,
    },
    tuple(['-ap', '--array-partition']): {
        'help': phelp_array_partition,
        'type': str,
        'default': DEFAULT_ARRAY_PARTITION,
    },
    tuple(['-aT', '--array-time']): {
        'help': phelp_array_time,
        'type': str,
        'default': DEFAULT_ARRAY_TIME,
    },
}

# ===========================================================================
# Program / serializer names + composed dicts
# ===========================================================================
L2D_VOTER_progname = 'L2D_VoterModel'
L2D_VOTER_progname_shrt = 'L2DVM'
L2D_VOTER_description = f"""
    Computational resources for the VoterModel dynamics on 2D
    lattices: {L2D_VOTER_progname}.py
"""
L2D_VOTER_srun_description = f"""Serializer for {L2D_VOTER_progname}.py"""

L2D_VOTER_args = {**L2D_args, **Voter_args}
L2D_VOTER_opt_args = {**Voter_opt_args, **Voter_save_args}
# Voter persists enabled observables itself (savedisk); the per-observable save
# gates ARE the saving options, so there is no Ising-style --save-results master
# switch or -sf checkpoint cadence.
L2D_VOTER_action_args = {**action_args_dict}
L2D_VOTER_srun_opt_args = {
    **Voter_srun_list_args,
    **srun_opt_args,
    **Voter_srun_array_args,
}
L2D_VOTER_srun_action_args = {**srun_action_args}
