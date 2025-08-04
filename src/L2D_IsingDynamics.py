from parsers.L2D_IsingDynamics import parse_arguments, parser
from kernels.L2D_IsingDynamics import *

def main():
    args = parse_arguments(parser)
    run_simulation(args)
    #
    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")

if __name__ == "__main__":
    main()
