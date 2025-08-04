from lrgsglib.core import *
from lrgsglib.config.plotlib import (
    imshow_colorbar_caxdivider,
    generate_maxpercdiff_colormap,
)
from lrgsglib.config.utils import move_to_rootf, width_interval, dv
from lrgsglib.config.nx_patches import (
    signed_spectral_layout,
    get_kth_order_neighbours,
)

#
move_to_rootf(print_tf=True)

import sys

sys.setrecursionlimit(2500)


side = 10
eigenmode = 0
#
theLattice = Lattice2D(  #
    side1=side,
    geometry="squared",
)
SLRG_obj = SignedLaplacianAnalysis(  #
    sg=theLattice,
    initCond="all_1",  # f'ground_state_{eigenmode}'
    pflip=0.15,
    t_steps=10,
    no_obs=200,
)
#
SLRG_obj.init_weights()
SLRG_obj.flip_random_fract_edges()
SLRG_obj.compute_k_eigvV()
SLRG_obj.find_ising_clusters()
SLRG_obj.mapping_nodes_to_clusters()

mVsT = []
f = open('data/p.txt', 'a+')
for T in tqdm(np.linspace(0.001, 3, num=20)):
    magn_i = []
    for i in range(200):
        print(f"temperature {T:.3g}, replica: {i}", end='\r')
        # termalizzazione
        _, _ = SLRG_obj.run_ising_dynamics(nstepsIsing=200, T=T)
        magn, _ = SLRG_obj.run_ising_dynamics(
            nstepsIsing=10, T=T, save_magnetization=True
        )
        avgm_icluster = np.mean(
            [
                np.mean(SLRG_obj.magn_array_save[t][SLRG_obj.biggestIsing_cl])
                for t in range(len(SLRG_obj.magn_array_save))
            ]
        )

        magn_i.append(avgm_icluster)
        # misurare la magnetizzazione dell'isola maggiore,
        # mediarla dopo la termalizzazione
        pass
    mVsT.append([T, np.mean(magn_i)])
    print([T, np.mean(magn_i)])
    f.write("{}\t{}\n".format(T, np.mean(magn_i)))
f.close()
    

    # mediare tu tutte le realizzazioni
# diagramma di fase m VS T
