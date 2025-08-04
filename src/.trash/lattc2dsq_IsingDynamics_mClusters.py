from lattc2dsq_IsingDynamics_mClusters_Parser import *
#
# per il programma C salva <m>, <m^2>
# ogni volta sovrascrivi
#
args = parser.parse_args()
#
import_on = False
if not args.only_graphs:
    grphNdir = f"{DEFAULT_GRAPH_OUTDIR}N={args.L**2:d}/"
    grphpth = f"{grphNdir}{args.graph_filename}_p={args.p:.3g}.pkl"
    if os.path.exists(grphpth):
        import_on = True
#
stdFnameSFFX = (
    "" if args.graph_filename == DEFAULT_GRAPH_NAME else args.graph_filename
)
#
sqlatt = Lattice2D(
    side1=args.L,
    geometry="squared",
    pflip=args.p,
    import_on=import_on,
    stdFnameSFFX=stdFnameSFFX,
)
#
if not import_on:
    sqlatt.flip_random_fract_edges()
    sqlatt.export_graph()
    sqlatt.export_adj_bin()

ising_dyn = IsingDynamics(
    sg=sqlatt,
    T=args.T,
    ic=args.initial_cond,
    runlang="C1",
    NoClust=args.Nclust,
)
ising_dyn.init_ising_dynamics(randstring_OPT = False)
#
if args.only_graphs:
    # print("EXIT!", sys.argv)
    exit()
#
pathcl = (
    lambda i: (
        f"{ising_dyn.sg.path_ising}cl{i}_{ising_dyn.sg.std_fname}.bin"
    )
)
pathcl_list = [pathcl(i) for i in range(ising_dyn.NoClust)]
cluster_idxs = [os.path.exists(pathcl(i)) for i in range(ising_dyn.NoClust)]
#
if not all(cluster_idxs):
    ising_dyn.find_ising_clusters()
    ising_dyn.export_ising_clust()
#
for _ in range(args.number_of_averages):
    ising_dyn.run(tqdm_on=False)
