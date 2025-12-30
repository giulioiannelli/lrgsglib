from parsers.MCG_SlaplSpect import parse_arguments, parser
from kernels.MCG_SlaplSpect import *


def main():
    args = parse_arguments(parser)
    perform_spectral_calculations(args)

    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")


if __name__ == "__main__":
    main()
