from __future__ import annotations

import argparse

from lrgsglib.config.progargs import (
    L2D_SlaplSpect_progName,
    L2D_SlaplSpect_progNameShrt,
    L2D_SlaplSpect_srun_description,
    L2D_SlaplSpect_srun_optional_args_dict,
    L2D_SlaplSpect_args,
    L2D_SlaplSpect_optional_args_dict,
    L2D_SlaplSpect_action_args_dict,
    srun_action_args,
)
from parsers.shared import CustomHelpAction


inner_optional_args = {
    **L2D_SlaplSpect_optional_args_dict,
    **L2D_SlaplSpect_action_args_dict,
}

parser = argparse.ArgumentParser(
    description=L2D_SlaplSpect_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
    add_help=False,
)

parser.add_argument(
    '-h', '--help',
    action=CustomHelpAction,
    inner_prog_name=L2D_SlaplSpect_progName,
    inner_positional_args=L2D_SlaplSpect_args,
    inner_optional_args=inner_optional_args,
    help='show this help message and exit',
)

# Only sweep + slurm flags are added here; everything else passes through to
# L2D_SlaplSpect.py via *unknown args in the runner.
for key, opts in {**L2D_SlaplSpect_srun_optional_args_dict, **srun_action_args}.items():
    parser.add_argument(*key, **opts)


__all__ = [
    "parser",
    "L2D_SlaplSpect_progName",
    "L2D_SlaplSpect_progNameShrt",
]
