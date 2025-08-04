from lattc2dsq_percClusters_Parser import *
#
args = parser.parse_args()
#
merged_dict = Counter()

lattice = Lattice2D(args.L, pflip=args.p)
filename = f'{lattice.path_lrgsg}p={args.p:.3g}_na={args.number_of_averages}_{args.out_suffix}.pkl'
if os.path.exists(filename):
    os.exit(0)

for avg in range(args.number_of_averages):
    lattice = Lattice2D(args.L, pflip=args.p)
    lattice.flip_random_fract_edges()
    dist_dict = lattice.cluster_distribution_list()
    merged_dict += Counter(dist_dict)
    try:
        filenameold = f'{lattice.path_lrgsg}p={args.p:.3g}_na={avg}_{args.out_suffix}.pkl'
        os.remove(filenameold)
    except OSError:
        pass
    filename = f'{lattice.path_lrgsg}p={args.p:.3g}_na={avg+1}_{args.out_suffix}.pkl'
    with open(filename, 'wb') as file:
        pickle.dump(merged_dict, file)