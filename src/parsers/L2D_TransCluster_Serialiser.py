from lrgsglib.core import *
description = """
    Serialiser for Lattice2D_TransCluster.py
"""
# Default values for the optional parameters
DEFAULT_NUMBER_AVERAGES = 10**4  # Assuming a default value for demonstration
DEFAULT_OUT_SUFFIX = ""
DEFAULT_MODE = 'slanzarv_ordParam'
DEFAULT_GEO = "squared"
DEFAULT_CELL = "rand"
DEFAULT_PRINT = False
DEFAULT_EXEC = "no"
DEFAULT_mMB = 2**10
DEFAULT_MMB = 2**14
DEFAULT_FLOAT_TYPE = "float64"
DEFAULT_PREW = 0.
# Helpers for argparse descriptions
HELP_print = f"""
    Option to print the output of the Serialiser.
"""
HELP_exc = f"""
    Option to exec the output of the Serialiser.
"""
HELP_mMB = f"""
    Minimum MB quantity to be allocated for the single process | default={DEFAULT_mMB}
"""
HELP_MMB = f"""
    Maximum MB quantity to be allocated for the single process | default={DEFAULT_MMB}
"""
HELP_prew = f"""
    Probability of rewiring edges in the lattice. | default={DEFAULT_PREW}
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
HELP_mode = f"""
    mode, phaseTr or pCluster | default='{DEFAULT_MODE}'
"""
HELP_o = f"""
    Suffix for the output file name (optional) | default='{DEFAULT_OUT_SUFFIX}'
"""
HELP_t = f"""
    Floating point arithmetic depth. | default='{DEFAULT_FLOAT_TYPE}'
"""
Lattice2D_TransCluster_progName = "L2D_TransCluster"
Lattice2D_TransCluster_progNameShrt = "L2D"
# Setup the argument parser
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Required parameters
parser.add_argument(
    "--print",
    help=HELP_print,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--exec",
    help=HELP_exc,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "-mMB", "--slanzarv_minMB",
    help=HELP_mMB,
    default=DEFAULT_mMB,
    type=int,
)
parser.add_argument(
    "-MMB", "--slanzarv_maxMB",
    help=HELP_MMB,
    default=DEFAULT_MMB,
    type=int,
)
parser.add_argument(
    "--prew",
    default=DEFAULT_PREW,
    help=HELP_prew,
    type=float,
)
parser.add_argument(
    "-g",
    "--geometry",
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
    "-nA",
    "--number_of_averages",
    default=DEFAULT_NUMBER_AVERAGES,
    help=HELP_nA,
    type=int,
)
parser.add_argument(
    "-o",
    "--out_suffix",
    default=DEFAULT_OUT_SUFFIX,
    help=HELP_o,
    type=str,
)
parser.add_argument(
    "--mode",
    default=DEFAULT_MODE,
    help=HELP_mode,
    type=str,
)
parser.add_argument(
    "-t",
    "--float_type",
    default=DEFAULT_FLOAT_TYPE,
    help=HELP_t,
    type=str,
)

