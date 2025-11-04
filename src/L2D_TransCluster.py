from parsers.L2D_TransCluster import *
#
args = parser.parse_args()
#
side = args.L
p = args.p
geo = args.geometry
cell = args.cell_type
mode = args.mode
navg = args.number_of_averages
sfreq = args.save_frequency if args.save_frequency else navg // 20
outsx = args.out_suffix
typf = args.float_type
prew = args.prew
outd = args.outdir
#
match typf:
    case 'float32':
        typf = np.float32
    case 'float64':
        typf = np.float64
    case _:
        raise ValueError("Invalid float type specified")
#
match mode:
    case 'pCluster':
        extout = PKL
    case 'ordParam':
        extout = TXT
    case _:
        raise ValueError("Invalid mode specified")
#
def get_geometry_func(cell: str):
    match cell:
        case 'rand' | 'randZERR' | 'randXERR':
            def geometry_func(lattice: Lattice2D):
                return lattice.nwDict[cell]['G']
        case _ if cell.startswith('ball'):
            radius = get_first_int_in_str(cell)
            def geometry_func(lattice: Lattice2D):
                return lattice.nwDict.get_links_rball(radius)
        case _:
            raise ValueError("Invalid cell specified")

    return geometry_func
#
def file_path_maker(mpath, mode=mode, ppath = p, napath = navg, spath = outsx,
                    ctpath = cell, extout = extout, prew = prew):
    match mode:
        case 'pCluster':
            extout = PKL
        case 'ordParam':
            extout = TXT
    prewStr = f"prew={prew:.3g}" if prew != 0. else ""
    strName = '_'.join(filter(None, [mode, f"p={ppath:.3g}", ctpath, 
                     f"na={napath}", prewStr, spath]))
    return os.path.join(mpath, strName) + extout
#
geometry_func = get_geometry_func(cell)
testLattice = Lattice2D(side, pflip=p, geo=geo, sgpathn=args.outdir)
mpath = {'pCluster': testLattice.path_lrgsg, 
         'ordParam': testLattice.path_phtra.name}
filename = file_path_maker(mpath[mode])
if os.path.exists(filename):
    exit(f"File {os.path.split(filename)[1]} already exists.")
#
#
#
match mode:
    case 'pCluster':
        merged_dict = Counter()
        #
        nAvgDone = 0
        try:
            fnameExists = glob.glob(f"{file_path_maker(mpath[mode], napath='', 
                                                spath='', extout='')}*")[0]
            merged_dict = pk.load(open(fnameExists, 'rb'))   
            if outsx:
                avgIdx = -2
            else:
                avgIdx = -1
            nAvgDone = os.path.splitext(fnameExists.split('_')[avgIdx])[0]
            nAvgDone = int(re.search(r'\d+', nAvgDone).group())
            fnameOld = fnameExists
        except:
            fnameOld = file_path_maker(mpath[mode], napath=0)
        nAvgNeed = navg - nAvgDone
        #
        for current_period in range((nAvgNeed // sfreq) + bool(nAvgNeed % sfreq)):
            batch_size = min(nAvgNeed - current_period * sfreq, sfreq)
            for _ in range(batch_size):
                l = Lattice2D(side, pflip=p, geo=geo, init_nw_dict=True)
                l.flip_sel_edges(geometry_func(l))
                #
                l.compute_k_eigvV(typf=typf)
                dist_dict = l.get_cluster_distribution()
                merged_dict += dist_dict
            navgCurr = nAvgDone + (current_period + 1) * sfreq
            fnameNew = file_path_maker(mpath[mode], napath=navgCurr)
            try:
                os.rename(fnameOld, fnameNew)
            except FileNotFoundError or OSError:
                pass
            with open(fnameNew, "wb") as f:
                pk.dump(merged_dict, f)
            fnameOld = fnameNew
    case 'ordParam':
        Pinf = []
        Pinf_var = []
        neglinks = 0
        avg1 = 0
        for cont, avg in enumerate(range(navg)):
            l = Lattice2D(side, 
                          pflip=p, 
                          geo=geo, with_positions=True, 
                          init_nw_dict=True, 
                          sgpathn=args.outdir)
            l.flip_sel_edges(geometry_func(l))
            #
            l.compute_k_eigvV(typf=typf)
            # try:
            #     l.compute_pinf()
            # except IndexError:
            #     continue
            #
            avg1 += 1
            neglinks += l.Ne_n
            #
            l.load_eigV_on_graph(binarize=True)
            l.make_clustersYN("eigV0", +1)
            if len(l.clustersY) > 1:
                smax2 = len(max(l.clustersY, key=len)) / (1.0 * l.N)
            else:
                smax2 = 1.0
            #             #
            Pinf.append(smax2)
            Pinf_var.append(smax2)
            data=[avg1,
                l.pflip,
                neglinks/avg1,
                np.mean(Pinf),
                np.mean(Pinf_var),
                np.std(Pinf)]
            #
            if (avg1 % sfreq == 0):
                try:
                    filenameold = file_path_maker(mpath[mode], napath=avg1-sfreq)
                    os.remove(filenameold)
                except OSError:
                    pass
                filename = file_path_maker(mpath[mode], napath=avg1)
                with open(filename, 'wb') as file:
                    np.savetxt(file, np.atleast_2d(data), fmt='%.7g')
        remains = navg % sfreq
        os.rename(file_path_maker(mpath[mode], napath=navg-remains),
                  file_path_maker(mpath[mode], napath=navg))