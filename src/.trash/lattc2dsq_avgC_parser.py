from lrgsglib.core import *
from lrgsglib.config.const import *
import warnings

# Your code here...

# Suppress RuntimeWarnings in a specific code block
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    # Code that generates RuntimeWarnings

# The rest of your code here...

#
description = """
    Compute average specific heat of a square lattice with random fraction of 
    flipped edges (-1) equal to p.
"""
#
DIR_GRAPH_DEFAULT = DEFAULT_DATA_LATTICE2DSQ_OUTDIR + "graphs/"
DEFAULT_GRAPH_NAME = DEFLIST_LATTICE2D_GEOABBRV[1]
DEFAULT_OUT_SUFFIX = ""
#
DEFAULT_nA = f"| default={DEFAULT_NUNMBER_AVERAGES}"
DEFAULT_gn = f"| default={DEFAULT_GRAPH_NAME}"
DEFAULT_o = f"| default={DEFAULT_OUT_SUFFIX}"
#
HELP_L = f"""
    (int) L of the square lattice
"""
HELP_p = f"""
    (float) fraction of edges to flip
"""
HELP_nA = f"""
    (int) number of averages to compute {DEFAULT_nA:->10}
"""
HELP_gn = f"""
    (str) file name of the graph container {DEFAULT_gn:->10}
"""
HELP_o = f"""
    (str) suffix to be given to the output {DEFAULT_o:->10}
"""
#
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
parser.add_argument(
    "-nA",
    "-navg",
    "--number_of_averages",
    default=DEFAULT_NUNMBER_AVERAGES,
    help=HELP_nA,
    type=int,
)
parser.add_argument(
    "-gn",
    "--graph_filename",
    default=DEFAULT_GRAPH_NAME,
    help=HELP_gn,
    type=str,
)
parser.add_argument(
    "-o",
    "--out_suffix",
    default=DEFAULT_OUT_SUFFIX,
    help=HELP_o,
    type=str,
)