from pathlib import Path
from typing import Union, Tuple
from numpy import integer as np_int
from numpy import pi as np_pi
from dotenv import load_dotenv
#
from .lrgsg_env import *
# types
ColorType = Union[
    Tuple[int, int, int],
    Tuple[float, float, float],
    Tuple[int, int, int, int],
    Tuple[float, float, float, float],
    str
]
GeneralInteger = (int, np_int)
# extensions
BIN = ".bin"
GML = ".gml"
LOG = ".log"
MP4 = ".mp4"
PDF = ".pdf"
PKL = ".pkl"
TXT = ".txt"
XML = ".xml"
# paths
load_dotenv()
#
PATHNLLIB: str = Path(LRGSG_LLIB).name
#
PATHDATA = Path(LRGSG_DATA)#
#
#
PATHNDATA = PATHDATA.name
#
PATHNCLDT = 'cluster_data'
PATHNGRPH = 'graph'
PATHNISNG = 'ising'
PATHNLRGS = 'lrgsg'
PATHNPAPR = 'paper'
PATHNPHTR = 'phtra'
PATHNPLOT = 'plot'
PATHNSPEC = 'spect'
PATHNVM = 'voter'
PATHNCP = 'cntct'
PATHNSRW = 'srw'
# Graph-related paths (created by SignedGraph)
PATHN_GRAPH_LIST = [PATHNGRPH, PATHNLRGS, PATHNPHTR, PATHNSPEC]
# Dynamics-related paths (created by respective dynamics classes)
PATHN_DYNAMICS_LIST = [PATHNISNG, PATHNVM, PATHNCP, PATHNSRW]
# All paths (for backward compatibility)
PATHN_LIST = PATHN_GRAPH_LIST + PATHN_DYNAMICS_LIST
#
PATHPLOT =  PATHDATA / Path(PATHNPLOT)
PATHPAPER = PATHPLOT / Path(PATHNPAPR)
PATHCLDAT = Path(PATHNDATA, PATHNCLDT)
#
# default values
LB_PFLIP = 0.
UB_PFLIP = 1.
#
BSP_VAL = [-1, +1]
BSP_STORE_MODE = ''
BSP_STORE_MODES = [BSP_STORE_MODE, 'persistent', 'sequential']
BSP_STORE_FREQ = 0
BSP_RUN_MODE_C = 'C'
BSP_RUN_MODE_PYTHON = 'Python'
BSP_RUN_MODE = BSP_RUN_MODE_PYTHON
BSP_RUN_MODES_C_LIST = [BSP_RUN_MODE_C, "c"]
BSP_RUN_MODES_C_DICT = {bsp: BSP_RUN_MODE_C for bsp in BSP_RUN_MODES_C_LIST}
BSP_RUN_MODES_PY_LIST = [BSP_RUN_MODE_PYTHON, "py", "Py", "PY", "Python"]
BSP_RUN_MODES_PY_DICT = {bsp: BSP_RUN_MODE_PYTHON 
                         for bsp in BSP_RUN_MODES_PY_LIST}
BSP_RUN_MODES = BSP_RUN_MODES_PY_LIST + BSP_RUN_MODES_C_LIST
BSP_RUN_MODES_DICT = {**BSP_RUN_MODES_PY_DICT, **BSP_RUN_MODES_C_DICT}
# Signed Graph default values
SG_PFLIP = 0.0
SG_IMPORT_ON = False
SG_IMPORT_EIGV = False
SG_INIT_NW_DICT = False
SG_INIT_WVAL = 1.0
SG_EXPORT_M = 'pk'
SG_LOAD_M = 'pk'
SG_LIST_REPR = ['H']
SG_REPR = 'G'
SG_GRAPHINT_REPR = SG_REPR
SG_GRAPHGEO_REPR = 'H'
SG_ERRMSG_MAXEIGVIDX = """The maximum eigenvalue index is out of bounds."""
SG_ERRMSG_NW_DICT = f"Inheriting class must have attribute 'nwContainer'"
SG_ERRMSG_NFLIP = """The probability of flipping an edge times the 
    number of edges is < 1, then no edges would be
    flipped. No flip will be performed."""
SG_ERRMSG_PFLIP = f""" pflip must be between {LB_PFLIP} and {UB_PFLIP},
    inclusive."""
SG_ERRMSG_NOCLUST_GTR_NCB = """Requested number of Cluster files is bigger than 
    the one in selected topology."""
SG_ERRMSG_ZEROCLUST = """No clusters were found."""
SG_WARNMSG_HASATTR = """The object does not has attribute."""
SG_WARNMSG_NOCLUST = """No clusters were found in the graph."""
# erdos renyi default values
WS_PHTABB = "ws"
WS_ONREP = 'G'
WS_STDFN = ""
WS_SGPATH = ""
# fully connected default values
FC_N = 100
FC_PHTABB = "fc"
FC_STDFN = ""
FC_SGPATH = ""
# scs generalized nn default values  
SCS_PHTABB = "scs_nn"
SCS_STDFN = ""
SCS_SGPATH = ""
# erdos renyi default values
ER_PHTABB = "er"
ER_ONREP = 'G'
ER_STDFN = ""
ER_SGPATH = ""
# Dorogovtsev-Goltsev-Mendes default values
DGM_PHTABB = "dgm"
DGM_ONREP = 'G'
DGM_STDFN = ""
DGM_SGPATH = ""
# 2D Lattice default values
L2D_BEND_POS = False
L2D_FBCV = 1.
L2D_PREW = 0.
L2D_PBC = True
L2D_STDFN = ""
L2D_SIDE1 = 32
L2D_SIDE2 = 0
L2D_SGPATH = ""
L2D_ONREP = 'G'
L2D_ONLY_CONST_MODE = False
L2D_PHTABB = 'l2d_'
L2D_GEO_TRI = 'triangular'
L2D_GEO_SQR = 'squared'
L2D_GEO_HEX = 'hexagonal'
L2D_GEO_SQRSW = 'squared_sw'
L2D_GEO_TRISW = 'triangular_sw'
L2D_GEO_SQOCT = 'octagonal_sqr'
L2D_GEO_KGM = 'kagome'
L2D_GEO_TRIHEX = 'tri_hexagonal'
L2D_GEO_TRI_SHRT = 'tri'
L2D_GEO_SQR_SHRT = 'sqr'
L2D_GEO_HEX_SHRT = 'hex'
L2D_GEO_SQRSW_SHRT = 'sqr_sw'
L2D_GEO_TRISW_SHRT = 'tri_sw'
L2D_GEO_SQOCT_SHRT = 'oct_sqr'
L2D_GEO_KGM_SHRT = 'kgm'
L2D_GEO_TRIHEX_SHRT = 'tri_hex'
L2D_P_C_LIST = [
    0.146,
    0.103,
    0.065,
    float("nan"),
    float("nan"),
    float("nan"),
    0.095,
    0.041,
] 
L2D_WITH_POS = False
# 
L2D_GEO = L2D_GEO_SQR
# 
L2D_GEO_LIST = [
    L2D_GEO_TRI,
    L2D_GEO_SQR,
    L2D_GEO_HEX,
    L2D_GEO_SQRSW,
    L2D_GEO_TRISW,
    L2D_GEO_SQOCT,
    L2D_GEO_KGM,
    L2D_GEO_TRIHEX,
]
L2D_GEO_SHRT_LIST = [
    L2D_GEO_TRI_SHRT,
    L2D_GEO_SQR_SHRT,
    L2D_GEO_HEX_SHRT,
    L2D_GEO_SQRSW_SHRT,
    L2D_GEO_TRISW_SHRT,
    L2D_GEO_SQOCT_SHRT,
    L2D_GEO_KGM_SHRT,
    L2D_GEO_TRIHEX_SHRT,
]
L2D_SINGLE_CELL_LIST = ['single', 'singleXERR', 'singleZERR']
L2D_RAND_CELL_LIST = ['rand', 'randXERR', 'randZERR']
#
L2D_P_C_DICT = {g: p for g,p in zip(L2D_GEO_LIST, L2D_P_C_LIST)}
L2D_PATH_DICT = {g: L2D_PHTABB + g for g in L2D_GEO_LIST}
L2D_ZIP_GEO_SHRT = list(zip(L2D_GEO_SHRT_LIST, L2D_GEO_LIST))
L2D_GEO_SHRT_DICT = {s: a for a,s in L2D_ZIP_GEO_SHRT}
L2D_SHRT_GEO_DICT = {a: s for a,s in L2D_ZIP_GEO_SHRT}
#
L2D_WARNMSG_GEO = """The selected geometry of the 2D lattice is not available. 
    Setting it to 'squared' for a 2d regular grid."""
L2D_ERRMSG_GEO = """Invalid side value for hexagonal lattice. In order to 
    implement PBC on hexagonal lattice you need to provide an even value 
    for side1 and side2."""
#
L3D_DIM0 = 16
L3D_DIM = tuple(L3D_DIM0 for _ in range(3))
L3D_PBC = True
L3D_FBCV = 1.
L3D_PDIL = 0.
L3D_SGPATH = ""
L3D_STDFN = ""
L3D_PHTABB = "l3d_"
L3D_GEO_SC = 'simple_cubic'
L3D_GEO_BCC = 'body_centered'
L3D_GEO_FCC = 'face_centered'
L3D_GEO_SC_SHRT1 = 'cubic'
L3D_GEO_SC_SHRT = 'sc'
L3D_GEO_BCC_SHRT = 'bcc'
L3D_GEO_FCC_SHRT = 'fcc'
L3D_ONREP = 'G'
L3D_ONLY_CONST_MODE = False
L3D_WITH_POS = False
L3D_THETA = np_pi/6
L3D_PHI = np_pi/6
# 
L3D_GEO = L3D_GEO_SC
#
L3D_GEO_LIST = [L3D_GEO_SC, 
                        L3D_GEO_BCC, 
                        L3D_GEO_FCC]
L3D_GEO_SHRT_LIST = [L3D_GEO_SC_SHRT, 
                        L3D_GEO_BCC_SHRT, 
                        L3D_GEO_FCC_SHRT]
L3D_GEO_DICT = {a: g for a,g in
                 zip(L3D_GEO_SHRT_LIST+L3D_GEO_LIST, L3D_GEO_LIST*2)}

# Multispectral graph default values
MSG_DEFAULT_TYPE = "multiplicative_cascade"
MSG_P1 = 0.8
MSG_P2 = 0.6
MSG_P3 = 0.6
MSG_P4 = 0.8
MSG_ITERATIONS = 7
MSG_FRACTION = 0.4
MSG_STDFN = ""
MSG_SGPATH = ""
MSG_ONLY_CONST_MODE = False
MSG_PHTABB = "msg_"
MSG_PATH = "msg_multiplicative_cascade"
MSG_PATHS = {MSG_DEFAULT_TYPE: MSG_PATH}

# MultiplicativeCascade constants
MC_STDFN = ""
MC_SGPATH = "msg_mc"
MC_PHTABB = "mc_"

# Vicsek constants
VCK_STDFN = ""
VCK_SGPATH = "msg_plv"
VCK_PHTABB = "vck_"

# DiracComb constants
DCOMB_STDFN = ""
DCOMB_SGPATH = "msg_dcomb"
DCOMB_PHTABB = "dcomb_"

# DiracBrush constants
DBRUSH_STDFN = ""
DBRUSH_SGPATH = "msg_dbrush"
DBRUSH_PHTABB = "dbrush_"
L3D_PATH_DICT = {a: L3D_PHTABB + g for a,g in 
                 zip(L3D_GEO_SHRT_LIST+L3D_GEO_LIST, L3D_GEO_LIST*2)}
#
L3D_ZIP_GEO_SHRT = list(zip(L3D_GEO_SHRT_LIST, L3D_GEO_LIST))
L3D_GEO_SHRT_DICT = {s: a for a,s in L3D_ZIP_GEO_SHRT}
L3D_SHRT_GEO_DICT = {a: s for a,s in L3D_ZIP_GEO_SHRT}
#
L3D_WARNMSG_GEO = """The selected geometry of the 3D lattice is not available. 
    Setting it to 'sc' for a 3d regular grid."""
#
DEFAULT_RECURSION_LIMIT = 1024**2
LRSG_ENTROPY_STEP = 1000
DEFAULT_ENTROPY_LEXPONENT = -3
DEFAULT_ENTROPY_HEXPONENT = 6
DEFAULT_ISING_NOCLUST = 5

# SignedRW walker defaults (see statsys/SignedRW/SignedWalker.py)
DEFAULT_SRW_N_WALKERS = 200
DEFAULT_SRW_SEED = 42
DEFAULT_SRW_COVERAGE_FRAC = 0.20
DEFAULT_SRW_MAX_N_CROSS = 14
DEFAULT_SRW_RULE = 'absorb'
DEFAULT_SRW_START = 'random'
DEFAULT_SRW_X_NODE = 'reflect'
SRW_RULES = ('absorb', 'kill', 'sticky')
SRW_START_PROTOCOLS = ('random', 'fixed', 'center')
SRW_X_NODE_BEHAVIORS = ('reflect', 'absorb')
SRW_OVERLAP_KINDS = ('l2_raw', 'l2_norm', 'l2_rescaled', 'fidelity', 'tv', 'jaccard')

DEFAULT_NUNMBER_AVERAGES = 100
DEFAULT_SPIKE_THRESHOLD = 0.05
DEFAULT_MAX_THRESHOLD = 2 * DEFAULT_SPIKE_THRESHOLD

DEFAULT_LOG_DIR = LRGSG_LOG
DEFAULT_P_FSTR_FMT = '.3g'

# Logging configuration constants
import logging as _logging
LOG_DEFAULT_LEVEL = _logging.WARNING
LOG_MAX_BYTES = 10 * 1024 * 1024  # 10 MB
LOG_BACKUP_COUNT = 3

DEFAULT_MAX_DIGITS_ROUND_SIGFIG = 18

COUNT_XERR_PATTERNS = True

# CHL-SF04 protein-TMD substrate defaults (see utils/recon/protein/substrate.py)
DEFAULT_ER_DEGREE = 10
DEFAULT_BA_M = 4
DEFAULT_WS_K = 6
DEFAULT_WS_P = 0.1
DEFAULT_RGG_RADIUS_3D = 0.18
DEFAULT_SBM_BLOCKS = 3
DEFAULT_SBM_PIN = 0.05
DEFAULT_SBM_POUT = 0.005

# CHL-SF04 sweep grids
PFLIP_SWEEP_GRID = (0.0, 0.05, 0.10, 0.20, 0.35, 0.50)
K_FRAC_GALLERY = (0.05, 0.20, 0.50)
K_FRAC_Q3 = (0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75)
N_SEEDS_PFLIP = 5
CHL_SF04_MASTER_SEED = 0xC4D14
CHL_SF04_RES_MIN = 40
CHL_SF04_RES_MAX = 200
CHL_SF04_CORPUS_SIZE = 10
Q3_EMERGENCE_THRESHOLD = 0.8
# Fixed feature dimensionality for all substrates in CHL-SF04 — every protein's
# coordinate vector is zero-padded to length N_FEATURES, every substrate is
# sized so its node count is >= N_FEATURES. The choice mirrors the original
# 32x32=1024 baseline (Iannelli & Villegas 2025 reproduces image recon at the
# same dimensionality) and makes per-substrate bases directly comparable.
CHL_SF04_N_FEATURES = 1024
CHL_SF04_LOG_MSE_FLOOR = 1e-4
