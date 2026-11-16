"""Topological percolation on random regular graphs (TransCluster v2).

    python src/RRG_TransCluster.py N k --p-linsp '(0, 0.5, 26)' [options]

``-tid t`` runs one tile (a SLURM array task); the default ``-tid -1`` runs
every tile in this process. Flow, run directory and exit codes:
kernels/RRG_TransCluster.py.
"""

import sys

from threadpoolctl import threadpool_limits

from kernels.RRG_TransCluster import resolve_run, run_transcluster_rrg
from lrgsglib.config.progargs.defs.generic import DEFAULT_SRUN_THREADS
from parsers.RRG_TransCluster import build_parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    command = sys.argv if argv is None else [parser.prog, *argv]
    try:
        run = resolve_run(args, command)
    except ValueError as err:
        parser.error(str(err))
    with threadpool_limits(limits=DEFAULT_SRUN_THREADS):
        return run_transcluster_rrg(run)


if __name__ == "__main__":
    sys.exit(main())
