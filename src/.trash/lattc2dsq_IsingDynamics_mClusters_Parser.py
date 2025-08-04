from lrgsglib.core import *
#
description = """
    Compute average magnetization of a square lattice with prescribed adjacency 
    matrix.
"""
#
DIR_GRAPH_DEFAULT = DEFAULT_DATA_LATTICE2DSQ_OUTDIR + "graphs/"
DEFAULT_GRAPH_NAME = DEFLIST_LATTICE2D_GEOABBRV[1]
DEFAULT_OUT_SUFFIX = ""
DEFAULT_INITCON = "ground_state_0"
#
DEFAULT_nA = f"| default={DEFAULT_NUNMBER_AVERAGES}"
DEFAULT_nC = f"| default={DEFAULT_ISING_NOCLUST}"
DEFAULT_gn = f"| default={DEFAULT_GRAPH_NAME}"
DEFAULT_o = f"| default={DEFAULT_OUT_SUFFIX}"
#
HELP_L = f"""
    (int) L of the square lattice
"""
HELP_p = f"""
    (float) fraction of edges to flip
"""
HELP_T = f"""
    (float) temperature of the system
"""
HELP_Nclust = f"""
    (int) the number of eigenstate whose biggest cluster has to be analyzed {DEFAULT_nC:->10}
"""
HELP_nA = f"""
    (float) temperature of the system {DEFAULT_nA:->10}
"""
HELP_gn = f"""
    (str) file name of the graph container {DEFAULT_gn:->10}
"""
HELP_o = f"""
    (str) file name of the graph container {DEFAULT_o:->10}
"""
#
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "L",
    help=HELP_L,
    type=int,
)
parser.add_argument(
    "T",
    help=HELP_T,
    type=float,
)
parser.add_argument(
    "p",
    help=HELP_p,
    type=float,
)
parser.add_argument(
    "-nC",
    "--Nclust",
    default=DEFAULT_ISING_NOCLUST,
    help=HELP_Nclust,
    type=int,
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
parser.add_argument(
    "-iC",
    "--initial_cond",
    default=DEFAULT_INITCON,
    type=str,
)
parser.add_argument(
    "--only_graphs",
    action=argparse.BooleanOptionalAction
)