from lrgsglib.core import *
plt.rcParams.update({'font.size': 30})

no_p = 10
no_rw = 5
pow2_m = 6
pow2_M = 14

lsN = 2**np.array([i for i in range(pow2_m, pow2_M, 2)])
lsL = np.sqrt(lsN)
lsp = {N: np.array([2 * N * p // (no_p) for p in range(1, no_p)]) for N in lsN}
lsNavg = 2**np.array([(pow2_M-i+1) for i in range(pow2_m, pow2_M, 2)])




for N, no_avg in zip(lsN, lsNavg):
    lsp = logpsace_prob_erconn(N)
    for p in lsp:
        for r in lsp:
            print(N, p, r)
            entropy_full = []
            Cspec_full = []
            fig, ax = plt.subplots(figsize=(12, 9))
            ax2 = ax.twinx()
            ax.set_ylabel(r'$1-S$')
            ax2.set_ylabel(r'$C\log N$')
            ax.set_xlabel(r'$\tau$')
            ax.set_xscale('log')
            for nvg in tqdm(range(no_avg)):
                G = nx.erdos_renyi_graph(N, p, seed=int(N*p), directed=False)
                if (no_avg-1 == nvg):
                    fignx, axnx = plt.subplots(figsize=(12, 9))
                    nx.draw(G, ax=axnx, node_size=300, with_labels = True)
                    fignx.savefig(f"{pltPath_Sm1C}ErdosReny_N={N}_p={p:.3g}_r={r:.3g}_eggraph{ePDF}", bbox_inches="tight")
                    plt.close('all')
                for e in G.edges():
                    G.add_edge(e[0], e[1], weight=1)
                    if np.random.random() < r:
                        G[e[0]][e[1]]['weight'] = -1
                [Sm1, dS1, VarL1, t1] = entropy(G, is_signed=True)
                dS1 /= np.max(dS1)
                t11 = (t1[1:] + t1[:-1]) / 2.
                entropy_full.append([t1, Sm1])
                Cspec_full.append([t11, dS1])
                pSm1, = ax.plot(t1, Sm1, '-')
                ax2.plot(t11, dS1, '--', c=pSm1.get_color())
            colors = ['r' if G[u][v]['weight'] == -1 else 'b' for u,v in G.edges()]
            fignx, axnx = plt.subplots(figsize=(12, 9))
            nx.draw(G, ax=axnx, edge_color=colors, node_size=300, with_labels = True)
            fignx.savefig(f"{pltPath_Sm1C}ErdosReny_N={N}_p={p:.3g}_r={r:.3g}_eggraph_r{ePDF}")
            plt.close(fignx)
            tSm1, avgdSm1 = averaged_Sm1(entropy_full)
            # tCs, avgdCspec = averaged_Sm1(Cspec_full)
            ax.plot(tSm1, avgdSm1, 'k-', lw=3)
            # ax.plot(tCs, avgdCspec, '--', c=pSm1.get_color())
            fig.savefig(f"{pltPath_Sm1C}Lattice2D/N={N}_p={p:.3g}_r={r:.3g}_navg={no_avg}{ePDF}")
            plt.close('all')