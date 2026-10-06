"""Command line of RRG_TransCluster.py (checks: kernels/RRG_TransCluster)."""

import argparse

from lrgsglib.config.progargs.RandomRegular import (
    RRG_args,
    RRG_TransCluster_description,
    RRG_TransCluster_grid_args,
    RRG_TransCluster_opt_args,
)


def build_parser() -> argparse.ArgumentParser:
    """Parser of RRG_TransCluster.py; exactly one control-grid flag."""
    parser = argparse.ArgumentParser(
        description=RRG_TransCluster_description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )
    for name, opts in RRG_args.items():
        parser.add_argument(name, **opts)
    grid = parser.add_mutually_exclusive_group(required=True)
    for flags, opts in RRG_TransCluster_grid_args.items():
        grid.add_argument(*flags, **opts)
    for flags, opts in RRG_TransCluster_opt_args.items():
        parser.add_argument(*flags, **opts)
    return parser
