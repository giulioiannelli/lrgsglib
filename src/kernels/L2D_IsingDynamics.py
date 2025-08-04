from lrgsglib import *
from .L2D import *
from .IsingDynamics import *

def run_simulation(args):
    ic_gs = args.init_cond.startswith('ground_state') or args.init_cond.startswith('gs')
    number = int(args.init_cond.split('_')[-1]) if ic_gs else 0
    val = ConditionalPartitioning(args.val)
    for _ in range(args.number_of_averages):
        lattice = Lattice2D(**initialize_l2d_dict_args(args))
        if lattice.init_nw_dict:
            lattice.flip_sel_edges(lattice.nwDict[args.cell_type]['G'])
        else:
            lattice.flip_random_fract_edges()
        lattice.compute_k_eigvV(number+1)
        lattice.load_eigV_on_graph(number, binarize=True)
        lattice.make_clustersYN(f'eigV{number}', val=val)
        #
        NoClust = lattice.handle_no_clust(NoClust=args.NoClust)
        isingDictArgs = initialize_ising_dict_args(args, get_out_suffix(args), NoClust)
        #
        if NoClust is None:
            continue
        #
        isdy = IsingDynamics(lattice, **isingDictArgs)
        isdy.init_ising_dynamics(exName=isdy.run_id)
        match args.runlang:
            case 'C4'|'C1':
                lattice.export_ising_clust(exName=isdy.run_id, 
                                        NoClust=NoClust)
                if args.runlang == 'C4':
                    lattice._export_eigV(number, exName=isdy.run_id, verbose=args.verbose)
        isdy.run(verbose=args.verbose)
        if args.remove_files:
            isdy.remove_run_c_files(remove_stderr=True)
            lattice.remove_exported_files()