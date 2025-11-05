# -*- coding: utf-8 -*-
# ##############################################################################
#
#   SCS_TransCluster.py
#   (Generalized SCS Network Transient Cluster Module)
#
#   Thin orchestration layer delegating parsing and heavy computations to
#   dedicated submodules, mirroring the structure of the lattice programs.
#
# ##############################################################################

"""Entry-point for the generalized SCS transient cluster analysis program."""

from parsers.SCS_TransCluster import *  # noqa: F401,F403
from kernels.SCS_TransCluster import run_transcluster


# ##############################################################################
#   Main Program Logic
# ##############################################################################

def main() -> None:
    """Parse CLI arguments and dispatch the transient cluster computation."""
    args = parse_arguments(parser)
    run_transcluster(args)


if __name__ == "__main__":
    main()
