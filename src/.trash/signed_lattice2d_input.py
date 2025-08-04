from lrgsglib.core import *
warnings.filterwarnings("ignore")
#
description = """
    Compute entropy of the density matrix pdf and the Signed Laplacian spectrum
    associated to a Random Bond Ising Model.
    Must be executed from parent folder "/path/to/LRG-Signed", i.e. called as
    python3 src/signed_lattice2dsq_input.py ARGUMENTS
"""
hPars_L = """
    (int) system linear size (L = sqrt(N))
"""
hPars_p = """
    (float) fraction of spins to flip 
"""
default_nA = "| default=1000"
hPars_nA = f"""
    (int) number of averages {default_nA:->10}
"""
default_st = "| default=1500"
hPars_st = f"""
    (int) number of steps in the specific heat sampling {default_st:->10}
"""
default_geom = "| default=squared"
hPars_geom = f"""
    (str) geometry of the 2d lattice {default_geom:->10}
"""
default_t1 = "| default=-2"
hPars_t1 = f"""
    (float) lower exponent of the time scale: {default_t1:->10}
"""
default_t2 = "| default=5"
hPars_t2 = f"""
    (float) lower exponent of the time scale: {default_t2:->10}
"""
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('L',
                    help=hPars_L,
                    type=int)
parser.add_argument('p',
                    help=hPars_p,
                    type=float)
parser.add_argument('-t1', '--lowert_exponent',
                    default=DEFAULT_ENTROPY_LEXPONENT,
                    help=hPars_t1,
                    type=float)
parser.add_argument('-t2', '--highert_exponent', 
                    default=DEFAULT_ENTROPY_HEXPONENT,
                    help=hPars_t2,
                    type=float)
parser.add_argument('-nA', '--number_of_averages',
                    default=DEFAULT_NUNMBER_AVERAGES,
                    help=hPars_nA,
                    type=int)
parser.add_argument('-st', '--sampling_steps',
                    # choices=range(1, int(1e10)),
                    default=DEFAULT_ENTROPY_STEPS,
                    help=hPars_st,
                    type=int)
parser.add_argument('-geo', '--geometry',
                  default='squared',
                #   choices=['squared', 'triangular', 'hexagonal'],
                  help=hPars_geom,
                  type=str)
args = parser.parse_args()
#
STEPS = args.sampling_steps
GEOMETRY = args.geometry
#
sideL = args.L
N = sideL * sideL
nA = args.number_of_averages
pflip = args.p
#
if not os.path.isdir((path := f"{datPath_l2d(GEOMETRY)}N={N}_navg={nA}/")):
    os.makedirs(path)
#
lattice = Lattice2D(side1 = sideL, geometry = GEOMETRY)
#
savename = lambda idstr: f"{path}{idstr}_p={pflip:{pflip_fmt}}{eBIN}"
if os.path.exists(savename('Sm1_avg')):
    exit()
Sm1AvgFile = open(savename('Sm1_avg'), "wb")
# slspecAvgFile = open(savename('slspec_avg'), "wb")
tau_maxFile = open(savename('tau_max'), "wb")
tau_maxFile0 = open(savename('tau_max0'), "wb")
Sm1lst = []
slspeclst = []
#
for nr in range(nA):
    random.seed(nr)
    while True:
        SLRG_obj = SignedLaplacianAnalysis(#
            sg = lattice,
            pflip = pflip,
            taumex = args.lowert_exponent,
            tauMex = args.highert_exponent,
            steps = STEPS
        )
        SLRG_obj.upd_graph_matrices()
        SLRG_obj.init_weights()
        SLRG_obj.flip_random_fract_edges()
        SLRG_obj.computeC()
        index = np.where(np.diff(SLRG_obj.Cspe) > DEFAULT_SPIKE_THRESHOLD)
        if not index[0].size:
            break
    SLRG_obj.compute_taumax_array()
    # slspeclst.append(SLRG_obj.slspectrum)
    Sm1lst.append(SLRG_obj.Sm1)
    tau_maxFile.write(bytes(SLRG_obj.taumax))
    tau_maxFile0.write(bytes(SLRG_obj.taumax0))
Sm1arr = np.array(Sm1lst).mean(axis=0)
# slspecarr = np.array(slspeclst).mean(axis=0)
#
Sm1AvgFile.write(bytes(Sm1arr))
# slspecAvgFile.write(bytes(slspecarr))
#
Sm1AvgFile.close()
# slspecAvgFile.close()
tau_maxFile.close()
tau_maxFile0.close()