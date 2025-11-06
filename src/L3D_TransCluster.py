"""Entry-point for the 3D transient cluster analysis program."""

from parsers.L3D_TransCluster import *  # noqa: F401,F403
from parsers.shared import parse_arguments

from kernels.L3D_TransCluster import run_transcluster


def main() -> None:
    """Parse CLI arguments and dispatch the transient cluster computation."""

    args = parse_arguments(parser)
    run_transcluster(args)


if __name__ == "__main__":
    main()
