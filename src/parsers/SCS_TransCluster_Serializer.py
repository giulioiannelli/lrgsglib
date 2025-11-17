from __future__ import annotations

import argparse
import sys

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


class CustomHelpAction(argparse.Action):
    """Custom help action that shows both serializer and inner program options."""
    
    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest, default=default, nargs=0, help=help)
    
    def __call__(self, parser, namespace, values, option_string=None):
        # Print serializer help
        parser.print_help()
        
        # Print inner program options
        print(f"\n{SCS_TransCluster_progName}.py options:")
        print(f"  Additional options passed to the inner {SCS_TransCluster_progName}.py program.")
        print(f"  These are not parsed by the serializer but passed through directly.\n")
        
        # Create a temporary parser for the inner program to display its help
        inner_parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            add_help=False,
        )
        
        # Add positional arguments
        for key, opts in SCSGeneralized_args.items():
            inner_parser.add_argument(key, **opts)
        
        # Add optional arguments
        inner_optional_args = {
            **SCSGeneralized_optional_args_dict,
            **SCS_TransCluster_optional_args_dict,
            **SCS_TransCluster_action_args_dict,
        }
        for key, opts in inner_optional_args.items():
            inner_parser.add_argument(*key, **opts)
        
        # Print the inner program's help (just the arguments section)
        inner_parser.print_help()
        
        sys.exit(0)


parser = argparse.ArgumentParser(
    description=SCS_TransCluster_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
    add_help=False,  # We'll add custom help
)

# Add custom help argument
parser.add_argument(
    '-h', '--help',
    action=CustomHelpAction,
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
