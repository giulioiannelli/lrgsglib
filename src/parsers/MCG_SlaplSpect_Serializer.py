from __future__ import annotations

import argparse

from lrgsglib.config.progargs import (
    MCG_SlaplSpect_srun_description,
    MCG_SlaplSpect_srun_optional_args_dict,
    MCG_MCGSS_progname,
    MCG_MCGSS_progname_shrt,
    srun_action_args,
    srun_opt_args,
    MC_args,
    MC_opt_args,
    MCG_SlaplSpect_optional_args_dict,
    MCG_SlaplSpect_action_args_dict,
)
from parsers.shared import CustomHelpAction

# Inner program arguments (what MCG_SlaplSpect.py accepts)
inner_optional_args = {
    **MC_opt_args,  # MC-specific options (--stochastic, --periodic, --variant)
    **MCG_SlaplSpect_optional_args_dict,  # Spectral options (--mode, --bins_count, etc.)
    **MCG_SlaplSpect_action_args_dict,
}

# Create parser
parser = argparse.ArgumentParser(
    description=MCG_SlaplSpect_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
    add_help=False,
)

# Add custom help argument to show serializer and inner program options
parser.add_argument(
    '-h', '--help',
    action=CustomHelpAction,
    inner_prog_name=MCG_MCGSS_progname,
    inner_positional_args=MC_args,
    inner_optional_args=inner_optional_args,
    help='show this help message and exit',
)

# Add serializer-specific argument group
serializer_group = parser.add_argument_group(
    'Serializer options',
    'Parameters for batch job submission and parameter sweeps',
)
for key, opts in {**MCG_SlaplSpect_srun_optional_args_dict, **srun_opt_args, **srun_action_args}.items():
    serializer_group.add_argument(*key, **opts)


__all__ = [
    "parser",
    "MCG_MCGSS_progname",
    "MCG_MCGSS_progname_shrt",
]
