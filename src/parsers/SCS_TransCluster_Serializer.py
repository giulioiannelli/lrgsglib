from __future__ import annotations

import argparse

from lrgsglib.config.progargs import (
    SCS_TransCluster_srun_description,
    SCS_TransCluster_srun_optional_args_dict,
    SCS_TransCluster_progName,
    SCS_TransCluster_progNameShrt,
    srun_action_args,
    SCSGeneralized_args,
    SCSGeneralized_optional_args_dict,
    SCS_TransCluster_optional_args_dict,
    SCS_TransCluster_action_args_dict,
)
from parsers.shared import CustomHelpAction


inner_optional_args = {
    **SCSGeneralized_optional_args_dict,
    **SCS_TransCluster_optional_args_dict,
    **SCS_TransCluster_action_args_dict,
}


parser = argparse.ArgumentParser(
    description=SCS_TransCluster_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
    add_help=False,  # We'll add custom help
)

# Add custom help argument that includes the inner program options
parser.add_argument(
    '-h', '--help',
    action=CustomHelpAction,
    inner_prog_name=SCS_TransCluster_progName,
    inner_positional_args=SCSGeneralized_args,
    inner_optional_args=inner_optional_args,
    help='show this help message and exit'
)

# Add serializer-specific arguments
serializer_group = parser.add_argument_group(
    'Serializer options',
    'Parameters for batch job submission and parameter sweeps'
)
for key, opts in {**SCS_TransCluster_srun_optional_args_dict, **srun_action_args}.items():
    serializer_group.add_argument(*key, **opts)


__all__ = [
    "parser",
    "SCS_TransCluster_progName",
    "SCS_TransCluster_progNameShrt",
]
