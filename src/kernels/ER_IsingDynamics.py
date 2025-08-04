from lrgsglib.core import *

MAX_ITER =1e5

def initialize_er_dict_args(args):
    return dict(n=args.N, p=args.p, sgpathn=args.workdir, pflip=args.pflip, init_nw_dict=True)

def initialize_ising_dict_args(args, out_suffix):
    return dict(T=args.T, ic=args.init_cond, runlang=args.runlang, NoClust=args.NoClust, rndStr=True, 
                out_suffix=out_suffix, id_string=args.in_suffix)

def get_out_suffix(args, ic_gs, number):
    return args.out_suffix or "gs"+str(number)+args.cell_type if ic_gs else args.init_cond

def run_simulation(args, erDictArgs, isingDictArgs, number, remove_files):
    count = 0
    while count < args.number_of_averages and count < MAX_ITER:
        er = ErdosRenyi(**erDictArgs)
        er.flip_sel_edges(er.nwDict[args.cell_type]['G'])
        er.compute_k_eigvV(k=number+1)
        er.load_eigV_on_graph(which=number, binarize=True)
        try:
            er.make_clustersYN(f'eigV{number}', -1)
        except:
            continue
        isdy = IsingDynamics(er, **isingDictArgs)
        isdy.init_ising_dynamics()
        er._export_edgel_bin(exName=isdy.run_id)
        er.export_ising_clust()
        isdy.run(verbose=False, thrmSTEP=args.thrmsteps)
        if remove_files:
            isdy.remove_run_c_files(remove_stderr=True)
            er.remove_edgl_file()
        count += 1
