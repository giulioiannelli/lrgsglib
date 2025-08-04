from lrgsglib.core import *

description = """
    This program computes the cluster distribution size of the sign fluctuation inside the signed Laplacian operator of a square lattice. It requires specifying the size of the lattice and the flipping probability. Optionally, it can average the results over multiple runs.
"""

# Default values for the optional parameters
DEFAULT_NUMBER_AVERAGES = 100  # Assuming a default value for demonstration
DEFAULT_OUT_SUFFIX = ""

# Helpers for argparse descriptions
HELP_L = """
    Size (L) of the square lattice.
"""

HELP_p = """
    Flipping probability (p) for edges.
"""

HELP_nA = f"""
    Number of averages to compute (optional) | default={DEFAULT_NUMBER_AVERAGES}
"""

HELP_o = f"""
    Suffix for the output file name (optional) | default='{DEFAULT_OUT_SUFFIX}'
"""

# Setup the argument parser
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required parameters
parser.add_argument(
    "L",
    help=HELP_L,
    type=int,
)

parser.add_argument(
    "p",
    help=HELP_p,
    type=float,
)

# Optional parameters
parser.add_argument(
    "-nA", "--number_of_averages",
    default=DEFAULT_NUMBER_AVERAGES,
    help=HELP_nA,
    type=int,
)

parser.add_argument(
    "-o", "--out_suffix",
    default=DEFAULT_OUT_SUFFIX,
    help=HELP_o,
    type=str,
)