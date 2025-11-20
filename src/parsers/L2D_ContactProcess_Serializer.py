from __future__ import annotations

import argparse

from lrgsglib.config.progargs import (
    L2D_CPROC_srun_description,
    L2D_CPROC_srun_optional_args_dict,
    L2D_CPROC_progname,
    L2D_CPROC_progname_shrt,
    srun_action_args,
    L2D_CPROC_args,
    L2D_CPROC_opt_args,
    L2D_CPROC_action_args,
)
from parsers.shared import CustomHelpAction


inner_optional_args = {
    **L2D_CPROC_opt_args,
    **L2D_CPROC_action_args,
}


parser = argparse.ArgumentParser(
    description=L2D_CPROC_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
    add_help=False,
)

parser.add_argument(
    "-h",
    "--help",
    action=CustomHelpAction,
    inner_prog_name=L2D_CPROC_progname,
    inner_positional_args=L2D_CPROC_args,
    inner_optional_args=inner_optional_args,
    help="show this help message and exit",
)

serializer_group = parser.add_argument_group(
    "Serializer options",
    "Parameters for batch job submission and parameter sweeps",
)

for key, opts in {**L2D_CPROC_srun_optional_args_dict, **srun_action_args}.items():
    serializer_group.add_argument(*key, **opts)


__all__ = [
    "parser",
    "L2D_CPROC_progname",
    "L2D_CPROC_progname_shrt",
]
