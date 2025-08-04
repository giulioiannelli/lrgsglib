from lrgsglib.core import *
description = """
    pCluster: Cluster distribution in signed Laplacian on lattices. 
    ordParam: Calculates order parameter and susceptibility.
"""
# Default values for the optional parameters
DEFAULT_NUMBER_AVERAGES = 1000  # Assuming a default value for demonstration
DEFAULT_SAVING_FREQUENCY = 0
DEFAULT_OUT_SUFFIX = ""
DEFAULT_MODE = "ordParam"
DEFAULT_GEO = "sc"
DEFAULT_CELL = "rand"
DEFAULT_FLOAT_TYPE = "float64"
DEFAULT_PDIL = 0.
DEFAULT_MU = 0.
DEFAULT_SIGMA = 0.
DEFAULT_EDGE_WEIGHT = 'flip'
# Helpers for argparse descriptions
HELP_L = """
    Size of the 3D lattice.
"""
HELP_p = """
    Flipping probability for edges.
"""
HELP_pdil = """
    Dilution probability of links. | default={DEFAULT_PDIL}
"""
HELP_geo = f"""
    Geometry of the lattice. | default={DEFAULT_GEO}
"""
HELP_cell = f"""
    Topological defect class: 'rand', 'randXERR', 
    'randZERR', 'ball_<R>' with type(<R>)=int. | default= '{DEFAULT_CELL}'
"""
HELP_nA = f"""
    Number of averages to compute | default={DEFAULT_NUMBER_AVERAGES}
"""
HELP_mode = f"""
    mode, ordParam or pCluster | default='{DEFAULT_MODE}'
"""
HELP_o = f"""
    Suffix for the output file name | default='{DEFAULT_OUT_SUFFIX}'
"""
HELP_svfreq = f"""
    Saving frequency for the data. By default 1/20 of the number of averages 
    | default='{DEFAULT_NUMBER_AVERAGES // 20}'
"""
HELP_t = f"""
    Floating point arithmetic depth. | default='{DEFAULT_FLOAT_TYPE}'
"""
HELP_mu = f"""
    Mean of the normal distribution for the edge weights. | default='{DEFAULT_MU}'
"""
HELP_sigma = f"""
    Standard deviation of the normal distribution for the edge weights. | default='{DEFAULT_SIGMA}'
"""
HELP_edge_weight = f"""
    Edge weight mode for the links. | default='{DEFAULT_EDGE_WEIGHT}'
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
    "--pdil",
    default=DEFAULT_PDIL,
    help=HELP_pdil,
    type=float
)
parser.add_argument(
    "--mu",
    default=DEFAULT_MU,
    help=HELP_mu,
    type=float
)
parser.add_argument(
    "--sigma",
    default=DEFAULT_SIGMA,
    help=HELP_sigma,
    type=float
)
parser.add_argument(
    "--edge_weight",
    default=DEFAULT_EDGE_WEIGHT,
    help=HELP_edge_weight,
    type=str
)
parser.add_argument(
    "--mode",
    default=DEFAULT_MODE,
    help=HELP_mode,
    type=str,
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
    "-n",
    "--number_of_averages",
    default=DEFAULT_NUMBER_AVERAGES,
    help=HELP_nA,
    type=int,
)
parser.add_argument(
    "--save_frequency",
    default=DEFAULT_SAVING_FREQUENCY,
    help=HELP_svfreq,
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
    "-t",
    "--float_type",
    default=DEFAULT_FLOAT_TYPE,
    help=HELP_t,
    type=str,
)

