from lrgsglib.core import *

ssize_list = [2**k for k in range(4, 8)]
nrep = [2**(17-2*k) for k in range(4, 8)]
nrep_range = [range(nr) for nr in nrep]
pval_list = [0.01, 0.05, 0.09, 0.098, 0.099, 0.100, 0.101, 0.102, 0.103, 0.104, 0.15, 0.2, 0.3]
print('list of sizes', ssize_list,
      '\nlist of flip probs', pval_list,
      '\nlist of replicas', nrep)

d_lmin = []

for i,L in enumerate(ssize_list):
    nedges = 2*L**2
    path = f"{datPath_lminl2d}N={L*L}_navg={nrep[i]}/"
    for p in pval_list:
        savename = f"{path}p={p:.3g}.txt"
        if os.path.exists(savename):
            continue
        if not os.path.isdir(path):
            os.makedirs(path)
        f = open(savename, "x")
        f.close()
        lmin = []
        for nr in tqdm(nrep_range[i], desc=f"replicas for L={L}, p={p}"):
            G = nx.grid_2d_graph(L, L, periodic=True)
            ransample = random.sample(range(nedges), int(p*nedges))
            #
            all_weights = {e: 1 for e in G.edges()}
            neg_weights = {e: -1 for i,e in enumerate(G.edges()) if i in ransample}
            #
            nx.set_edge_attributes(G, values=all_weights, name='weight')
            nx.set_edge_attributes(G, values=neg_weights, name='weight')
            #
            slapl = get_graph_lapl(G)
            try:
                lmin_tmp, _ = eigsh(slapl, k=1, which='SM', tol=1e-10)
                lmin.append(lmin_tmp)
            except ArpackNoConvergence:
                pass
        lminarr = np.array(lmin).flatten()
        d_lmin.append({'L': L, 'p':p, 'lmin': lminarr})
        
        np.savetxt(savename, lminarr)