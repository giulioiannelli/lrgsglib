"""SLURM arrays of RRG_TransCluster.py tiles (TransCluster v2).

    python src/RRG_TransCluster_Serializer.py --N-list 1000 4000 --k-list 3 \\
        --p-linsp '(0, 0.5, 26)' -na 200 [--print | --exec | --status |
        --collect] [--only-missing] [--verify]

Options it does not know go to every task unchanged. Planning, scripts and
actions: kernels/RRG_TransCluster_Serializer.py.
"""

import sys

from kernels.RRG_TransCluster_Serializer import (
    ACTIONS,
    check_options,
    plan_runs,
)
from parsers.RRG_TransCluster_Serializer import build_parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args, extras = parser.parse_known_args(argv)
    command = sys.argv if argv is None else [parser.prog, *argv]
    try:
        check_options(args, parser)
        plans = plan_runs(args.N_list, args.k_list, extras)
    except ValueError as err:
        parser.error(str(err))
    return ACTIONS[args.action](plans, args, command)


if __name__ == "__main__":
    sys.exit(main())
