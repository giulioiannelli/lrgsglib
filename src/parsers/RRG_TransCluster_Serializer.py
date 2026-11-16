"""Command line of RRG_TransCluster_Serializer.py.

Only the serializer's own options; every other option is forwarded to the
tasks (``parse_known_args``). ``-h`` lists the program's options too.
"""

import argparse

from lrgsglib.config.progargs.defs.TransCluster import DEFAULT_TC_SRUN_ACTION
from lrgsglib.config.progargs.RandomRegular import (
    RRG_TC_TILE_FLAGS,
    RRG_TransCluster_grid_args,
    RRG_TransCluster_opt_args,
    RRG_TransCluster_progName,
    RRG_TransCluster_srun_actions,
    RRG_TransCluster_srun_args,
    RRG_TransCluster_srun_description,
)
from parsers.shared import CustomHelpAction


def build_parser() -> argparse.ArgumentParser:
    """Parser of RRG_TransCluster_Serializer.py; one action at most."""
    parser = argparse.ArgumentParser(
        description=RRG_TransCluster_srun_description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
        add_help=False,
    )
    forwarded = {**RRG_TransCluster_grid_args, **RRG_TransCluster_opt_args}
    parser.add_argument(
        "-h",
        "--help",
        action=CustomHelpAction,
        inner_prog_name=RRG_TransCluster_progName,
        inner_positional_args={},  # N and k come from --N-list / --k-list
        inner_optional_args={
            flags: opts
            for flags, opts in forwarded.items()
            if flags != RRG_TC_TILE_FLAGS
        },
        help="show this help message and exit",
    )
    for flags, opts in RRG_TransCluster_srun_args.items():
        parser.add_argument(*flags, **opts)
    parser.set_defaults(action=DEFAULT_TC_SRUN_ACTION)
    actions = parser.add_mutually_exclusive_group()
    for flags, opts in RRG_TransCluster_srun_actions.items():
        actions.add_argument(*flags, **opts)
    return parser
