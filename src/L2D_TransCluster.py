# -*- coding: utf-8 -*-
# ##############################################################################
#
#   L2D_TransCluster.py
#   (2D Lattice Transient Cluster Module)
#
#   Part of the lrgsglib package for network topology analysis in lattice systems.
#
# ##############################################################################
#
#   Description:
#   ------------
#   Thin orchestration layer delegating the heavy lifting to the parser and
#   kernel submodules. The goal is to keep the entry-point clean and easily
#   readable, mirroring the structure adopted in ``L2D_Recon.py``.
#
# ##############################################################################

"""Entry-point for the 2D transient cluster analysis program."""

from parsers.L2D_TransCluster import *  # noqa: F401,F403
from parsers.shared import parse_arguments
from kernels.L2D_TransCluster import run_transcluster


# ##############################################################################
#   Main Program Logic
# ##############################################################################


def main() -> None:
    """Parse CLI arguments and dispatch the transient cluster computation."""
    args = parse_arguments(parser)
    run_transcluster(args)


if __name__ == "__main__":
    main()
