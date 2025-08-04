from lrgsglib.core import *

description = """ 
    Compute the order parameter of the transition and the susceptibility for a lattice setting.
"""

# Default values for the optional parameters
DEFAULT_NUMBER_AVERAGES = 1000 # Assuming a default value for demonstration
DEFAULT_OUT_SUFFIX = ""
DEFAULT_GEO = "squared"
DEFAULT_CELL = "single"

# Helpers for argparse descriptions
HELP_L = """
    Size (L) of the square lattice.
"""

HELP_p = """
    Flipping probability (p) for edges.
"""

HELP_geo = f"""
    Geometry of the lattice. | default={DEFAULT_GEO}
"""

HELP_cell = f"""
    Topological class of the defect: possible values 'single', 'square', 'triangle', 'hexagon', 'cross'. | default= '{DEFAULT_CELL}'
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
    "-g", "--geometry",
    default=DEFAULT_GEO,
    help=HELP_geo,
    type=str,
)

parser.add_argument(
    "-c", "--cell_type",
    default=DEFAULT_CELL,
    help=HELP_cell,
    type=str,
)

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