from __future__ import annotations

import argparse

from lrgsglib.config.progargs import (
    L3D_TransCluster_srun_description,
    L3D_TransCluster_srun_optional_args_dict,
    L3D_TransCluster_progName,
    L3D_TransCluster_progNameShrt,
    srun_action_args,
    L3D_TransCluster_args,
    L3D_TransCluster_optional_args_dict,
    L3D_TransCluster_action_args_dict,
)
from parsers.shared import CustomHelpAction


inner_optional_args = {
    **L3D_TransCluster_optional_args_dict,
    **L3D_TransCluster_action_args_dict,
}

parser = argparse.ArgumentParser(
    description=L3D_TransCluster_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
    add_help=False,
)

# Add custom help argument to show serializer and inner program options
parser.add_argument(
    '-h', '--help',
    action=CustomHelpAction,
    inner_prog_name=L3D_TransCluster_progName,
    inner_positional_args=L3D_TransCluster_args,
    inner_optional_args=inner_optional_args,
    help='show this help message and exit',
)

# NOTE: We only add serialiser-specific arguments here (list/linspace parameters
# and slurm options). All other arguments (geometry, cell_type, mode, etc.) should
# be passed directly to L3D_TransCluster.py via *unknown args, where they will be
# properly parsed by the underlying program's parser.

for key, opts in {**L3D_TransCluster_srun_optional_args_dict, **srun_action_args}.items():
    parser.add_argument(*key, **opts)


__all__ = [
    "parser",
    "L3D_TransCluster_progName",
    "L3D_TransCluster_progNameShrt",
]
