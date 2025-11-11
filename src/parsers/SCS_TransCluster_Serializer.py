from __future__ import annotations

import argparse

from lrgsglib.config.progargs import (
    SCS_TransCluster_srun_description,
    SCS_TransCluster_srun_optional_args_dict,
    SCS_TransCluster_progName,
    SCS_TransCluster_progNameShrt,
    srun_action_args,
)

parser = argparse.ArgumentParser(
    description=SCS_TransCluster_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
)

# NOTE: We only add serialiser-specific arguments here (list/linspace parameters
# and slurm options). All other arguments (workdir, backend, mode, etc.) should
# be passed directly to SCS_TransCluster.py via *unknown args, where they will be
# properly parsed by the underlying program's parser.

for key, opts in {**SCS_TransCluster_srun_optional_args_dict, **srun_action_args}.items():
    parser.add_argument(*key, **opts)


__all__ = [
    "parser",
    "SCS_TransCluster_progName",
    "SCS_TransCluster_progNameShrt",
]
