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
#
if not os.path.isdir((path := f"{datPath_l2d(GEOMETRY)}N={N}_navg={nA}/")):
    os.makedirs(path)
#
lattice = Lattice2D(side1 = sideL, geometry = GEOMETRY)
lattice.lsp_selection(lattice.default_dict_lsp(num_at=7))

for pflip in lattice.lsp:
    savename = lambda idstr, pflip: f"{path}{idstr}_p={pflip:{pflip_fmt}}{eBIN}"   
    print(pflip)
    mean_fract = []
    for nr in range(nA):
        random.seed(nr)
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
        eigv, eigV = eigsh(SLRG_obj.slp.astype(np.float64), k=1, which='SM')
        N = len(eigV)
        Nm = (eigV > 0).sum()
        mean_fract.append(np.abs(N-2*Nm)/N)
    np.savetxt(savename('fract', pflip), np.array(mean_fract))