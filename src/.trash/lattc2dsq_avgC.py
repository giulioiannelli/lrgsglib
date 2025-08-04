from lattc2dsq_avgC_parser import *
#
args = parser.parse_args()
#
stdFnameSFFX = (
    "" if args.graph_filename == DEFAULT_GRAPH_NAME else args.graph_filename
)
#
thresh = 2e-2
lst = []
_ = 0
while _ < args.number_of_averages:
    sqlatt = Lattice2D(
        side1=args.L,
        geometry="squared",
        pflip=args.p,
        stdFnameSFFX=stdFnameSFFX,
    )
    #
    sqlatt.flip_random_fract_edges()
    sla = SignedLaplacianAnalysis(sqlatt)
    sla.computeC()
    avgCval = medfilt(sla.Cspe)
    if all(avgCval[np.where(np.abs(np.diff(avgCval)) == avgCval[:-1])] < thresh):
        lst.append(avgCval)
        _ += 1
#     print(f'\r{_}', end="")
# print('\n')
file_path = f"{sqlatt.path_lrgsg}avgC_{args.number_of_averages}_p={sqlatt.pflip:.3g}_{args.out_suffix}.npz"
np.savez_compressed(file_path, *lst)