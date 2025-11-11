from __future__ import annotations

import argparse

from lrgsglib.config.progargs import (
    L3D_TransCluster_srun_description,
    L3D_TransCluster_srun_optional_args_dict,
    L3D_TransCluster_progName,
    L3D_TransCluster_progNameShrt,
    srun_action_args,
)

parser = argparse.ArgumentParser(
    description=L3D_TransCluster_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
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
