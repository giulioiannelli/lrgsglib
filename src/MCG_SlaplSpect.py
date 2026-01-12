import sys

# CRITICAL: Force unbuffered output for slurm logs
sys.stdout = open(sys.stdout.fileno(), 'w', buffering=1)  # Line buffered
sys.stderr = open(sys.stderr.fileno(), 'w', buffering=1)

print("[INIT] MCG_SlaplSpect starting...", flush=True)

from parsers.MCG_SlaplSpect import parse_arguments, parser
from kernels.MCG_SlaplSpect import *

print("[INIT] Imports complete", flush=True)


def main():
    args = parse_arguments(parser)
    print(f"[INIT] Args parsed: it={args.iterations}, fraction={args.fraction}, backend={args.backend}", flush=True)
    perform_spectral_calculations(args)

    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")


if __name__ == "__main__":
    main()
