from argparse import SUPPRESS

from .common import *
from .defs.RandomRegular import *
from .defs.TransCluster import *
from .phelp.RandomRegular import *
from .phelp.TransCluster import phelp_edge_weight, phelp_mu
from .SignedGraph import SG_opt_args

# ===========================================================================
# RRG_TransCluster.py: topological percolation on random regular graphs
# ===========================================================================
RRG_TransCluster_progName = "RRG_TransCluster"
RRG_TransCluster_progNameShrt = "RRGTC"
RRG_TransCluster_description = f"""
    Topological percolation of the signed-Laplacian ground state on random
    regular graphs, one tile of (rep, control value) units per call:
    {RRG_TransCluster_progName}.py
"""
RRG_TC_TILE_FLAGS = ("-tid", "--tile-id")  # the serializer sets it per task
RRG_args = {
    "N": {"help": phelp_rrg_N, "type": int},
    "k": {"help": phelp_rrg_k, "type": int},
}
# Control grid: exactly one of these (a required mutually exclusive group).
RRG_TransCluster_grid_args = {
    tuple(["--p-list"]): {
        "help": phelp_tc_p_list,
        "type": float,
        "nargs": "+",
    },
    tuple(["--p-linsp"]): {
        "help": phelp_tc_p_linsp,
        "type": parse_multiple_linspace,
    },
    tuple(["--sigma-list"]): {
        "help": phelp_tc_sigma_list,
        "type": float,
        "nargs": "+",
    },
    tuple(["--sigma-linsp"]): {
        "help": phelp_tc_sigma_linsp,
        "type": parse_multiple_linspace,
    },
}
RRG_TransCluster_opt_args = {
    # disorder
    tuple(["-c", "--cell_type"]): {
        **SG_opt_args[tuple(["-c", "--cell_type"])],
        "help": phelp_tc_cell,
        "choices": list(TC_CELL_FAMILY),
        "default": DEFAULT_TC_CELL,
    },
    tuple(["-ew", "--edge_weight"]): {
        "help": phelp_edge_weight,
        "type": str,
        "choices": list(TC_EDGE_WEIGHTS),
        "default": DEFAULT_TC_EDGE_WEIGHT,
    },
    tuple(["--mu"]): {
        "help": phelp_mu,
        "type": float,
        "default": DEFAULT_TC_MU,
    },
    tuple(["--ensemble"]): {
        "help": phelp_tc_ensemble,
        "type": str,
        "choices": list(TC_ENSEMBLES),
        "default": DEFAULT_TC_ENSEMBLE,
    },
    # statistics and tiling
    tuple(["-na", "--number_of_averages"]): {
        **SG_opt_args[tuple(["-na", "--number_of_averages"])],
        "default": DEFAULT_RRG_NAVG,
    },
    tuple(["--rep-start"]): {
        "help": phelp_tc_rep_start,
        "type": int,
        "default": DEFAULT_TC_REP_START,
    },
    tuple(["--ctrl-chunk"]): {
        "help": phelp_tc_ctrl_chunk,
        "type": int,
        "default": DEFAULT_TC_CTRL_CHUNK,
    },
    tuple(["--rep-chunk"]): {
        "help": phelp_tc_rep_chunk,
        "type": int,
        "default": DEFAULT_TC_REP_CHUNK,
    },
    RRG_TC_TILE_FLAGS: {
        "help": phelp_tc_tile_id,
        "type": int,
        "default": DEFAULT_TC_TILE_ID,
    },
    tuple(["--seed"]): {
        "help": phelp_rrg_seed,
        "type": int,
        "default": DEFAULT_RRG_SEED,
    },
    # topology and solver
    tuple(["-ge", "--graph_engine"]): {
        **graph_engine_opt_args[tuple(["-ge", "--graph_engine"])],
        "help": phelp_rrg_engine,
        "default": DEFAULT_RRG_ENGINE,
    },
    tuple(["--sampler"]): {
        "help": phelp_rrg_sampler,
        "type": str,
        "choices": list(RRG_SAMPLERS),
        "default": DEFAULT_RRG_SAMPLER,
    },
    tuple(["--eig-tol"]): {
        "help": f"{phelp_tc_eig_tol}; at most {TC_EIG_TOL_MAX:g}",
        "type": float,
        "default": DEFAULT_TC_EIG_TOL,
    },
    tuple(["--gap-every"]): {
        "help": phelp_tc_gap_every,
        "type": int,
        "default": DEFAULT_TC_GAP_EVERY,
    },
    # output
    tuple(["-wd", "--workdir"]): {
        **SG_opt_args[tuple(["-wd", "--workdir"])],
        "default": DEFAULT_RRG_WORKDIR,
    },
    tuple(["-os", "--out_suffix"]): {
        **common_opt_args[tuple(["-os", "--out_suffix"])],
        "help": phelp_rrg_tag,
    },
    tuple(["-v", "--verbose"]): action_args_dict[tuple(["-v", "--verbose"])],
}
# ===========================================================================
# RRG_TransCluster_Serializer.py: SLURM arrays of RRG_TransCluster tiles
# ===========================================================================
RRG_TransCluster_srun_description = f"""
    SLURM arrays of {RRG_TransCluster_progName}.py tiles, one run per (N, k).
    Options not listed here are forwarded to every task unchanged (see the
    {RRG_TransCluster_progName}.py options below); -tid is set per task.
"""
RRG_TransCluster_srun_args = {
    tuple(["--N-list"]): {
        "help": phelp_rrg_N_list,
        "type": int,
        "nargs": "+",
        "required": True,
    },
    tuple(["--k-list"]): {
        "help": phelp_rrg_k_list,
        "type": int,
        "nargs": "+",
        "required": True,
    },
    tuple(["--only-missing"]): {
        "help": phelp_tc_only_missing,
        "action": "store_true",
    },
    tuple(["--verify"]): {
        "help": phelp_tc_verify,
        "action": "store_true",
    },
    # resources
    tuple(["-ac", "--array-concurrent"]): {
        "help": phelp_tc_array_concurrent,
        "type": int,
        "default": DEFAULT_TC_ARRAY_CONCURRENT,
    },
    tuple(["-ap", "--array-partition"]): {
        "help": phelp_tc_array_partition,
        "type": str,
        "default": DEFAULT_TC_ARRAY_PARTITION,
    },
    tuple(["-aT", "--array-time"]): {
        "help": phelp_tc_array_time,
        "type": str,
        "default": DEFAULT_TC_ARRAY_TIME,
    },
    tuple(["--array-max"]): {
        "help": phelp_tc_array_max,
        "type": int,
        "default": DEFAULT_ARRAY_MAX_SIZE,
    },
    tuple(["--mem-mb"]): {
        "help": phelp_tc_mem_mb,
        "type": int,
        "default": DEFAULT_TC_MEM_MB,
    },
    # --collect
    tuple(["--boot"]): {
        "help": phelp_tc_boot,
        "type": int,
        "default": DEFAULT_TC_BOOT,
    },
    tuple(["--boot-seed"]): {
        "help": phelp_tc_boot_seed,
        "type": int,
        "default": DEFAULT_TC_BOOT_SEED,
    },
    tuple(["--peak-window"]): {
        "help": phelp_tc_peak_window,
        "type": float,
        "nargs": 2,
        "metavar": ("LO", "HI"),
    },
    tuple(["--legacy-txt"]): {
        "help": phelp_tc_legacy_txt,
        "action": "store_true",
    },
}
# Mutually exclusive modes; all store their name in args.action. The parser
# sets the default (DEFAULT_TC_SRUN_ACTION) once, so -h shows none per flag.
RRG_TransCluster_srun_actions = {
    tuple([f"--{name}"]): {
        "help": phelp_tc_srun_actions[name],
        "action": "store_const",
        "const": name,
        "dest": "action",
        "default": SUPPRESS,
    }
    for name in TC_SRUN_ACTIONS
}
