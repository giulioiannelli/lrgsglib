"""CLI entry point for the signed random-walker dynamics on 2D lattices."""

from parsers.L2D_SignedRW import parse_arguments, parser
from kernels.L2D_SignedRW import run_simulation
from lrgsglib.utils.tools.chronometer import Chronometer


def main() -> None:
    args = parse_arguments(parser)
    run_simulation(args)
    if getattr(args, "print_chrono", False):
        Chronometer.print_all_chronometers()
    if getattr(args, "verbose", False):
        print("Done!")


if __name__ == "__main__":
    main()
