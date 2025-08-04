from lrgsglib.core import *
#
path_suppinf = PATHPLOT / Path(PATHNPAPR, 'suppinf')
path_extras = path_suppinf / Path('extras')
path_isingXERR = path_extras / Path('isingXERR')
#
seed = 0xFADE
#
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("side", type=int, help="side of the sqlattice")  # This enforces float type
parser.add_argument("p", type=float, help="fraction to flip")  # This enforces float type
parser.add_argument("-T", type=float, help="temperature", default=0.5)  # This enforces float type
parser.add_argument("--cell", type=str, default='randXERR', help="cell type")

args = parser.parse_args()
#
side = args.side
pflip = args.p
T = args.T
cell = args.cell
#
ic_list = ["uniform"] + ["_".join(["gs", str(_)]) for _ in [0, 1, 2, 3, 4, 10, 50]]
#
l2d = Lattice2D(side, pflip=pflip, init_nw_dict=True, 
                    path_data=path_isingXERR, seed=seed)
l2d.flip_sel_edges(l2d.nwDict[cell]['G'])
l2d.compute_k_eigvV(50)
#
for ic in ic_list:
    fname = lambda who: Path(f"{who}_p={pflip:.3g}_T={T:.3g}_{ic}_{cell}.bin")
    pname_ene = l2d.path_ising / fname("ene")
    pname_sout = l2d.path_ising / fname("sout")
    if not os.path.exists(pname_ene) or not os.path.exists(pname_sout):
        isdy = IsingDynamics(l2d, T=T, ic=ic, runlang="C3", thrmSTEP=80, out_suffix=f"{ic}_{cell}")
        isdy.init_ising_dynamics()
        isdy.run()