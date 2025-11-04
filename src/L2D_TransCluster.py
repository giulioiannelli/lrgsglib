# -*- coding: utf-8 -*-
"""Entry point for the L2D Transient Cluster analysis program."""

from parsers.L2D_TransCluster import parse_arguments, parser
from kernels.L2D_TransCluster import execute_transcluster
from lrgsglib.utils.tools.Chronometer import Chronometer


def main() -> None:
    """Parse CLI arguments and execute the transient cluster workflow."""

    args = parse_arguments(parser)
    execute_transcluster(args)

    if getattr(args, "print_chrono", False):
        Chronometer.print_all_chronometers()
    if getattr(args, "verbose", False):
        print("Done!")


if __name__ == "__main__":
    main()
