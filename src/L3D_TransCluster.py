from parsers.L3D_TransCluster import *

def parse_args():
    args = parser.parse_args()
    return args

def get_file_path(mpath, mode, ppath, napath, spath, ctpath, extout, pdil, mu, sigma, edge_weight):
    pdilStr = f"pdil={pdil:.3g}" if pdil else None
    pflipStr = f"p={ppath:.3g}" if edge_weight == 'flip' else f"mu={mu:.3g}_sigma={sigma:.3g}"
    navgStr = f"na={napath}"
    listStr = [mode, pdilStr, pflipStr, ctpath, navgStr, spath]
    strName = '_'.join(filter(None, listStr))
    return os.path.join(mpath, strName) + extout

def get_geometry_func(cell):
    if cell in {'rand', 'randZERR', 'randXERR'}:
        return lambda lattice: lattice.nwDict[cell]['G']
    elif cell.startswith('ball'):
        radius = get_first_int_in_str(cell)
        return lambda lattice: lattice.nwDict.get_links_rball(radius)
    else:
        raise ValueError("Invalid cell specified")

def load_existing_data(filename):
    try:
        fnameExists = glob.glob(f"{filename}*")[0]
        merged_dict = pk.load(open(fnameExists, 'rb'))
        avgIdx = -2 if outsx else -1
        nAvgDone = int(re.search(r'\d+', os.path.splitext(fnameExists.split('_')[avgIdx])[0]).group())
        return merged_dict, nAvgDone, fnameExists
    except:
        return Counter(), 0, filename

def save_data(filename, data, mode):
    if mode == 'pCluster':
        with open(filename, "wb") as f:
            pk.dump(data, f)
    elif mode == 'ordParam':
        with open(filename, 'wb') as file:
            np.savetxt(file, np.atleast_2d(data), fmt='%.7g')

def process_pCluster(lattice, geometry_func, nAvgNeed, period, mpath, mode, fnameOld):
    merged_dict = Counter()
    for cp in range((nAvgNeed // period) + bool(nAvgNeed % period)):
        batch_size = min(nAvgNeed - cp * period, period)
        for _ in range(batch_size):
            l = Lattice3D(side, pflip=p, geo=geo, pdil=pdil, init_nw_dict=True)
            l.flip_sel_edges(geometry_func(l))
            merged_dict += l.get_cluster_distribution()
        navgCurr = nAvgDone + (cp + 1) * period
        fnameNew = get_file_path(mpath, mode, p, navgCurr, outsx, cell, extout, pdil, mu, sigma, edge_weight)
        try:
            os.rename(fnameOld, fnameNew)
        except (FileNotFoundError, OSError):
            pass
        save_data(fnameNew, merged_dict, mode)
        fnameOld = fnameNew

def process_ordParam(lattice, geometry_func, navg, sfreq, mpath, mode):
    Pinf, Pinf2, neglinks = [], [], 0
    for avg in range(navg):
        l = Lattice3D(side, pflip=p, geo=geo, pdil=pdil, init_nw_dict=True)
        if edge_weight == 'flip':
            l.flip_sel_edges(geometry_func(l))
        else:
            l.set_edges_random_normal(mu, sigma)
        neglinks += l.Ne_n
        l.compute_k_eigvV(typf=typf)
        try:
            l.compute_pinf()
        except IndexError:
            continue
        Pinf.append(l.Pinf)
        Pinf2.append(l.Pinf**2)
        if (avg + 1) % sfreq == 0:
            Pinf_tot = sum(Pinf) / (avg + 1)
            Pinf2_tot = sum(Pinf2) / (avg + 1)
            data = [avg + 1, l.pflip, neglinks / (avg + 1), Pinf_tot, Pinf2_tot, Pinf2_tot - Pinf_tot**2]
            try:
                os.remove(get_file_path(mpath, mode, p, avg + 1 - sfreq, '', cell, extout, pdil, mu, sigma, edge_weight))
            except OSError:
                pass
            save_data(get_file_path(mpath, mode, p, avg + 1, outsx, cell, extout, pdil, mu, sigma, edge_weight), data, mode)

args = parse_args()
side, p, pdil, mu, sigma, edge_weight, geo, cell, mode, navg, sfreq, outsx, typf = (
    args.L, args.p, args.pdil, args.mu, args.sigma, args.edge_weight, args.geometry, args.cell_type, args.mode, args.number_of_averages, args.save_frequency if args.save_frequency else args.number_of_averages // 20, args.out_suffix, args.float_type
)

if mode == 'pCluster':
    extout = '.pkl'
elif mode == 'ordParam':
    extout = '.txt'
else:
    raise ValueError("Invalid mode specified")
geometry_func = get_geometry_func(cell)
lattice = Lattice3D(dim=tuple(side for _ in range(3)), pflip=p, geo=geo)
mpath = {'pCluster': lattice.path_lrgsg, 'ordParam': lattice.path_phtra.name}[mode]
filename = get_file_path(mpath, mode, p, navg, outsx, cell, extout, pdil, mu, sigma, edge_weight)

if os.path.exists(filename):
    exit(f"File {os.path.split(filename)[1]} already exists.")
else:
    merged_dict, nAvgDone, fnameOld = load_existing_data(get_file_path(mpath, mode, p, '', '', cell, extout, pdil, mu, sigma, edge_weight))
    nAvgNeed = navg - nAvgDone
    lattice = Lattice3D(dim=tuple(side for _ in range(3)), pflip=p, geo=geo)
    if mode == 'pCluster':
        process_pCluster(lattice, geometry_func, nAvgNeed, sfreq, mpath, mode, fnameOld)
    elif mode == 'ordParam':
        process_ordParam(lattice, geometry_func, navg, sfreq, mpath, mode)
