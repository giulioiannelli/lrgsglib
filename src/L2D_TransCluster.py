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
            # Pinf.append(l.Pinf)
            # Pinf_var.append(l.Pinf_var)
            # l.make_clustersYN("eigV0", +1)
            # # print(l.clustersY)
            # # print(list(map(len, l.clustersY)))
            # nodelist_clust = list(sorted(map(len, l.clustersY), reverse=True))
            #
            l.load_eigV_on_graph(binarize=True)
            l.make_clustersYN("eigV0", +1)
            if len(l.clustersY) > 1:
                smax2 = len(max(l.clustersY, key=len)) / (1.0 * l.N)
            else:
                smax2 = 1.0
            # largest_component_subG = max(nx.connected_components(subG), key=len)
            # largest_component_clustersY = max(l.clustersY, key=len)
            # if set(largest_component_subG) != set(largest_component_clustersY) or avg == 1:
            #     print("Discrepancy found!")
            #     print(smax, smax2)
            #     print("pablo:", sorted(list(largest_component_subG)))
            #     print("giulio:", sorted(list(largest_component_clustersY)))
            #     plt.figure(figsize=(18, 5))
            #     # Plot 1: Network with labels and eigenstate overlay
            #     plt.subplot(1, 4, 1)
            #     node_colors = [l.G.nodes[node].get("eigV0", 0) for node in l.G.nodes]  # Use eigV0 attribute if it exists, otherwise 0
            #     nx.draw(l.G, pos=l.get_node_attributes(), node_color=node_colors, cmap=plt.cm.viridis, with_labels=True)
            #     plt.title("Network with Node Labels and Eigenstate Overlay")

            #     # Plot 2: Network colored with attribute `eigV0` of each node
            #     plt.subplot(1, 4, 2)
            #     giant_component_nodes = set(largest_component_clustersY)  # The largest connected component from `giulio`
            #     node_colors_giant = ["red" if node in giant_component_nodes else "blue" for node in l.G.nodes]
            #     nx.draw(l.G, pos=l.get_node_attributes(), node_color=node_colors_giant, cmap=plt.cm.viridis, with_labels=False)
            #     plt.title("Network Colored by Attribute `eigV0`")

            #     # Plot 3: Network with subgraph colored differently
            #     plt.subplot(1, 4, 3)
            #     subG = l.get_subgraph_from_nodes(largest_component_subG)  # Assuming this extracts the required subgraph
            #     subgraph_nodes = set(subG.nodes())
            #     full_node_colors = ["red" if node in subgraph_nodes else "blue" for node in l.G.nodes]
            #     nx.draw(l.G, pos=l.get_node_attributes(), node_color=full_node_colors, with_labels=False)
            #     plt.title("Network with Subgraph in Red and Rest in Blue")

            #     # Plot 2: Network colored with attribute `eigV0` of each node
            #     plt.subplot(1, 4, 4)

            #     nx.draw(l.G, pos=l.get_node_attributes(), node_color=list(l.get_node_attributes('eigV0').values()), cmap=plt.cm.viridis, with_labels=False)
            #     plt.title("Network Colored by Attribute `eigV0`")

            #     plt.tight_layout()
            #     plt.show()
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