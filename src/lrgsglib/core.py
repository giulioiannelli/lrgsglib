from .config import *
from .utils import *
from .nx_patches import *
from .statsys import *
#
# __all__ = [
#     "np",
#     "nx",
#     "plt",
#     "powerlaw",
#     "SignedLaplacianAnalysis",
#     "IsingDynamics",
#     "FullyConnected",
#     "Lattice2D",
# ]
#
# renormalization group for heterogenous network functions
class SignedLaplacianAnalysis:
    Sm1 = None
    VarL = None
    Cspe = None
    taumax = None
    taumax0 = None
    frames_dynsys = []
    Ising_clusters = []
    eigv = None
    ACCERR_LAPL_DYN = 1e-10
    MAXVAL_LAPL_DYN = 200
    #
    def __init__(
        self,
        sg: SignedGraph,
        taustep: int = LRSG_ENTROPY_STEP,
        taulex: float = -2,
        tauhex: float = 5,
        maxThresh: float = DEFAULT_MAX_THRESHOLD,
        nreplica: int = 0,
        initCond: str = "gauss_1",
        t_steps=10,
        no_obs=1,
        initspect: bool = True
    ) -> None:
        self.sg = sg
        self.initspect = initspect
        self.nreplica = nreplica
        self.taustep = taustep
        self.taulex = taulex
        self.tauhex = tauhex
        self.tTsS = np.logspace(self.taulex, self.tauhex, self.taustep)
        self.tTsC = .5 * (self.tTsS[1:] + self.tTsS[:-1])
        self.maxThresh = maxThresh
        self.initCond = initCond
        #
        self.t_steps = t_steps
        self.no_obs = no_obs
        if not hasattr(self.sg, "eigv") and self.initspect:
            self.__initSpectrum__()
    #
    def __initSpectrum__(self):
        self.sg.compute_k_eigvV(with_routine="numpy")
    #
    def computeS(self) -> None:
        w =  self.sg.eigv
        S = np.zeros(len(self.tTsS))
        VarL = np.zeros(len(self.tTsS))
        for i, tau in enumerate(self.tTsS):
            rhoTr = np.exp(-tau * w)
            Tr = np.nansum(rhoTr)
            rho = np.divide(rhoTr, Tr)
            avgrho = np.nansum(np.multiply(w, rhoTr)) / Tr
            av2rho = np.nansum(np.multiply(np.multiply(w, w), rhoTr)) / Tr
            S[i] = -np.nansum(rho * np.log(rho)) / np.log(self.sg.N)
            VarL[i] = av2rho - avgrho**2
        self.Sm1 = 1 - S
        self.VarL = VarL
    #
    def computeC(self) -> None:
        if self.Sm1 is None:
            self.computeS()
        self.Cspe = np.log( self.sg.N) * dv(self.Sm1, np.log(self.tTsS))

    #
    def compute_taumax_array(self, filter_gauss1d: int = 3) -> None:
        if self.Cspe is None:
            self.computeC()
        if filter_gauss1d:
            # if you dont do this sometimes the function is not enough sensible to detect the maxima
            Cspe = gaussian_filter1d(self.Cspe, sigma=filter_gauss1d) 
        else:
            Cspe = self.Cspe
        maxIdx = argrelextrema(Cspe, np.greater)[0]
        maxIdxCondition = self.Cspe[maxIdx] > self.maxThresh
        self.taumax = self.tTsC[maxIdx[maxIdxCondition]]
        self.taumax0 = self.tTsC[maxIdx[maxIdxCondition][0]]

    #
    def laplacian_dynamics_init(
        self,
        t_stepsMultiplier: int = 1,
        window_size=0,
        window_shift_x=0,
        window_shift_y=0,
        win_val=1,
    ):
        N =  self.sg.N
        #
        self.Deltat = 1.0 / self.t_steps
        self.simulationTime = N * t_stepsMultiplier
        self.sampling = np.logspace(
            0, np.log10(self.simulationTime), num=self.no_obs, dtype=int
        )
        #
        if not "ground_state" in self.initCond:
            self.sg.compute_k_eigvV()
        if self.initCond.startswith("ground_state"):
            self.eigenModeInit = int(self.initCond.split("_")[-1])
            self.sg.compute_k_eigvV(k=self.eigenModeInit + 1)
            self.field =  self.sg.eigV.T[self.eigenModeInit]
        elif self.initCond == "uniform_1":
            self.field = np.random.uniform(-1, 1, N)
        elif self.initCond == "delta_1":
            self.field = np.zeros(N)
            self.field[N // 2] = 1
        elif self.initCond == "gauss_1":
            self.field = np.random.normal(-1, 1, N)
        elif self.initCond.startswith("all"):
            self.initVal = float(self.initCond.split("_")[-1])
            self.field = self.initVal * np.ones(N)
        elif self.initCond.startswith("window"):
            s22 =  self.sg.side2 // 2
            wndwS = s22 - 1 if window_size > (s22 - 1) else window_size
            hS, hE = s22 - wndwS - 1, s22 + wndwS + 1
            initStatus = np.zeros(N)
            #
            if self.initCond.startswith("window_multiple"):
                self.nsquares = int(self.initCond.split("_")[-1])
                wndwSa = np.array([wndwS, -wndwS])
                sqTmp = np.random.randint(
                    wndwS,  self.sg.side1 - wndwS, size=(self.nsquares, 2)
                )
                result = np.column_stack((sqTmp, sqTmp + wndwS))
                result[:, [1, 2]] = result[:, [2, 1]]
                result = result.reshape(-1, 2, 2)
                sqIdx = np.concatenate(
                    [
                        np.concatenate(
                            [
                                [
                                    j + i *  self.sg.side1
                                    for j in range(*iSq[0])
                                ]
                                for i in range(*iSq[1])
                            ]
                        )
                        for iSq in result
                    ]
                )
                if type(win_val) is str:
                    if win_val == "gauss":
                        win_val = np.random.normal(0, 1, size=len(sqIdx))
                    elif win_val == "uniform":
                        win_val = np.random.uniform(-1, 1, len(sqIdx))
                initStatus[sqIdx] = np.ones(len(sqIdx)) * win_val
            else:
                shiftsX = [hS + window_shift_x, hE + window_shift_x]
                shiftsY = [hS + window_shift_y, hE + window_shift_y]
                sqIdx = np.concatenate(
                    [
                        [j + i *  self.sg.side1 for j in range(*shiftsX)]
                        for i in range(*shiftsY)
                    ]
                )
                if self.initCond == "window_1":
                    initStatus[sqIdx] = np.ones(len(sqIdx))
                elif self.initCond == "window_val":
                    initStatus[sqIdx] = np.ones(len(sqIdx)) * win_val
                elif self.initCond == "window_gauss_1":
                    initStatus[sqIdx] = np.random.normal(-1, 1, len(sqIdx))
                elif self.initCond.startswith("window_ground_state"):
                    self.eigenModeInit = int(self.initCond.split("_")[-1])
                    self.sg.compute_k_eigvV(k=self.eigenModeInit + 1)
                    initStatus = self.eigV.T[self.eigenModeInit]
                    outIdx = np.setxor1d(np.arange(N), sqIdx)
                    initStatus[outIdx] = np.random.uniform(
                        initStatus.min(), initStatus.max(), len(outIdx)
                    )
                else:
                    print("Error, no mode for init laplacian dynamic chosen.")
            self.field = initStatus
        #
        if  self.sg.pbc is False:
            L = int(np.sqrt(N))
            self.fixed_border_idxs = np.array(
                sorted(
                    [i for i in range(L)]
                    + [(L - 1) * L + i for i in range(L)]
                    + [i * L for i in range(1, L - 1)]
                    + [(i + 1) * L - 1 for i in range(1, L - 1)]
                )
            )
            self.field[self.fixed_border_idxs] =  self.sg.fbc_val

    #
    def run_laplacian_dynamics(self, rescaled=False, saveFrames=False):
        self.frames_dynsys = []
        #
        try:
            x = self.field
        except AttributeError:
            self.laplacian_dynamics_init()
            x = self.field
        #
        def stop_conditions_lapdyn(self, x_tm1, xx):
            ERRTOL = self.ACCERR_LAPL_DYN * np.ones( self.sg.N)
            C1 = (np.abs(x_tm1 / x_tm1.max() - xx / xx.max()) < ERRTOL).all()
            C2 = np.abs(np.log10(np.max(np.abs(xx)))) > self.MAXVAL_LAPL_DYN
            return C1, C2
        #
        if rescaled:
            if rescaled == "dynamic":
                lap = (
                    lambda t: np.exp(- self.sg.eigv[0] * t) *  self.sg.slp
                )
            else:
                self.sg.make_rescaled_signed_laplacian(rescaled)
                lap = lambda _:  self.sg.resLp
        else:
            lap = lambda _:  self.sg.slp
        #
        if not self.sg.pbc:
            def set_bc():
                x[self.fixed_border_idxs] = self.sg.fbc_val
        else:
            def set_bc():
                pass
        #
        if saveFrames:
            def save_frames(self, x, t):
                if t in self.sampling:
                    self.frames_dynsys.append(x)
        else:
            def save_frames(*_):
                pass
        #
        print("Beginning Laplacian dynamics.")
        #
        for t in tqdm(range(1, self.simulationTime)):
            save_frames(self, x, t)
            x_tm1 = x
            x = x - self.Deltat * lap(t) @ x  # + np.sqrt(Deltat)*noise
            set_bc()
            C1, C2 = stop_conditions_lapdyn(self, x_tm1, x)
            if C2 or C1:
                if C1:
                    print("Convergence reached.")
                if C2:
                    print("Max val. reached.")
                break
        self.field = x

    #
    def rescaled_field_regularization(self):
        status = self.field.reshape( self.sg.side1,  self.sg.side2)
        restatus = np.log10(np.max(status) - status)
        nnans = restatus[(restatus != np.inf) & (restatus != -np.inf)]
        self.restatus = np.nan_to_num(
            restatus, posinf=np.max(nnans), neginf=np.min(nnans)
        )

    #
    def make_animation_fromFrames(self, savename="output", fps=10, dpi=200):
        import matplotlib.animation as animation
        no_frames = len(self.frames_dynsys)
        print("# of frames: ", no_frames)
        #
        fig, ax = plt.subplots()
        #
        animate =  self.sg.make_animation(fig, ax, self.frames_dynsys)
        ani = animation.FuncAnimation(fig, animate, frames=no_frames)
        #
        fig.tight_layout()
        ani.save(
            f"{savename}{MP4}",
            writer=animation.FFMpegWriter(fps=fps),
            dpi=dpi,
        )
        plt.close(fig)


#                           ,((((((((((((((((((((((((.
#                      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                   @@@@@@                              @@@@@@
#                  @@@@@@                               *@@@@@/
#    /@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@,
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                #    TRASH CODE
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@,
#    @@@@@@@@@@@@@@    @@@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    &@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@&
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    .@@@@@@@@@@@@@
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    (@@@@@@@@@@@@@
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    @@@@@@@@@@@@@@
#     @@@@@@@@@@@@@.    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    @@@@@@@@@@@@@@
#     @@@@@@@@@@@@@/    @@@@@@@@@@@@@@     @@@@@@@@@@@@@&    @@@@@@@@@@@@@@
#     *@@@@@@@@@@@@@    @@@@@@@@@@@@@@     @@@@@@@@@@@@@(    @@@@@@@@@@@@@.
#      @@@@@@@@@@@@@    &@@@@@@@@@@@@@     @@@@@@@@@@@@@.    @@@@@@@@@@@@@
#      @@@@@@@@@@@@@    *@@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      &@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@#
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#       *@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@@    @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&
#          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&

#
# def lsp_read_values(folder_path, fpattern="Sm1_avg_p", sort=True):
#     file_pattern = rf"{fpattern}=(\d+\.\d+)"
#     value_pattern = r"p=(\d+\.\d+)"
#     # Get all files in the folder
#     files = os.listdir(folder_path)
#     # Filter files based on the pattern
#     file_names = [
#         file_name for file_name in files if re.match(file_pattern, file_name)
#     ]
#     # Extract values from file names
#     values = []
#     for file_name in file_names:
#         match = re.search(value_pattern, file_name)
#         if match:
#             value = float(match.group(1))
#             values.append(value)
#     # Sort the values if needed
#     if sort:
#         values.sort()
#     return np.array(values)

# if len(self.IsingIC) > len("ground_state"):
#     IICnoise = float(self.IsingIC.split("_")[-1])
#     if IICnoise < 0 or IICnoise > 1:
#         print("error")
#     self.s = self.s * np.random.choice(
#         [-1, 1], size=self.system.N, p=[IICnoise, 1 - IICnoise]
#     )

# self.distr = np.unique(result_array, return_counts=True)

#
# def find_ising_clusters(self):
#     if self.Ising_clusters:
#         print("exit function")
#         return
#     lnodes = list(self.system.H.nodes())
#     lnodes_tmp = lnodes[:]
#     #
#     self.system.compute_k_eigvV()
#     eigVbin = self.system.bin_eigV()
#     #
#     def recursive_search(seed, magn_i, clustertmp):
#         neighs = get_kth_order_neighbours(self.system.H, seed, 1)
#         neighs = np.array([e for e in neighs if e not in set(clustertmp)])
#         if not neighs.size:
#             return
#         samecluster = np.array(eigVbin[neighs] == magn_i)
#         if not samecluster.any():
#             return
#         neighs_samecluster = list(neighs[samecluster])
#         clustertmp.extend(neighs_samecluster)
#         for ss in neighs_samecluster:
#             recursive_search(ss, magn_i, clustertmp)
#     #
#     for i in lnodes:
#         if i not in lnodes_tmp:
#             continue
#         if not lnodes_tmp:
#             break
#         #
#         clustertmp = []
#         clustertmp.extend([i])
#         #
#         recursive_search(i, eigVbin[i], clustertmp)
#         lnodes_tmp = [e for e in lnodes_tmp if e not in set(clustertmp)]
#         self.Ising_clusters.append(clustertmp)
#     self.numIsing_cl = len(self.Ising_clusters)
#     self.Ising_clusters.sort(key=len, reverse=True)


# def run_ising_dynamics(
#     self,
#     T=0.1,
#     nstepsIsing=100,
#     IsingIC="uniform",
#     save_magnetization=False,
#     tqdm_on=False,
#     MODE='C'
# ):

#     self.init_ising_dynamics(IsingIC)
#     self.id_string_isingdyn = randstring()
#     if MODE == 'C':
#         output_file = open(f"src/lrgsglib/tmp_stuff/m_{self.id_string_isingdyn}.bin", 'wb')
#         self.s.tofile(output_file)
#         return

#     magn = []
#     ene = []
#     if save_magnetization:
#         self.magn_array_save = []
#     ising = self.IsingDynamics(self.system.N)
#     #
#     def boltzmann_factor(energy, temp):
#         return np.exp(-energy / temp)
#     #questo ci puo calcolare una sola volta
#     def neigh_weight_magn(m: NDArray, node_dict: dict) -> list:
#         return [w["weight"] * m[nn] for nn, w in node_dict.items()]

#     def neigh_ene(neigh: list) -> float:
#         return - np.sum(neigh) / neigh.__len__()

#     def calc_full_energy(m: NDArray, graph):
#         """Energy of a given configuration"""
#         return np.array(
#             [
#                 m[i] * neigh_ene(neigh_weight_magn(m, dict(graph.H[i])))
#                 for i in range(graph.N)
#             ]
#         ).sum()

#     # import time

#     # neigh_time1 = 0
#     # neigh_time2 = 0
#     # fliptime = 0

#     def flip_spin(node, m, graph, T):
#         # nonlocal neigh_time1
#         # nonlocal neigh_time2
#         # nonlocal fliptime
#         # nt1 = time.time()
#         m_flp = -m[node]
#         neigh = neigh_weight_magn(m, dict(graph.H[node]))
#         # ent1 = time.time()
#         # nt2 = time.time()
#         neighEn = neigh_ene(neigh)
#         E_old = m[node] * neighEn
#         E_new = m_flp * neighEn
#         # ent2 = time.time()
#         # et1 = time.time()
#         DeltaE = E_new - E_old
#         if DeltaE < 0:
#             m[node] = m_flp
#         elif np.random.uniform() < boltzmann_factor(DeltaE, T):
#             m[node] = m_flp
#         # net1 = time.time()
#         # neigh_time1 += ent1 - nt1
#         # neigh_time2 += ent2 - nt2
#         # fliptime += net1 - et1

#     flip_spin_all = np.vectorize(flip_spin)
#     flip_spin_all.excluded.add(1)
#     flip_spin_all.excluded.add(2)
#     flip_spin_all.excluded.add(3)
#     #
#     # s1 = time.time()
#     # # a = np.random.randint(0, self.system.N, nstepsIsing*self.system.N)
#     # t2_tot = 0.0
#     # t3_tot = 0.0

#     iterator = tqdm(range(nstepsIsing)) if tqdm_on else range(nstepsIsing)

#     for _ in range(nstepsIsing):
#         magn.append(np.sum(m))
#         ene.append(calc_full_energy(m, self.system))
#         # s2 = time.time()
#         sample = rd.sample(self.system.H.nodes(), self.system.N)
#         # e2 = time.time()
#         # t2_tot += e2 - s2
#         # s3 = time.time()
#         # for i in range(self.system.N):
#         #     # node =   # a[t*nstepsIsing + i]#sample[i]
#         #     flip_spin(sample[i], m, self.system, T)
#         flip_spin_all(sample, m, self.system, T)
#         # e3 = time.time()
#         # t3_tot += e3 - s3
#         if save_magnetization:
#             self.magn_array_save.append(np.array(m))
#     self.magn_array = m
#     # e1 = time.time()
#     # print("time for full:", e1 - s1)
#     # print("time for sampling:", t2_tot)
#     # print("time for flipping:", t3_tot)
#     # print("time for flipping_neig1:", neigh_time1)
#     # print("time for flipping_neig2:", neigh_time2)
#     # print("time for flipping_flip:", fliptime)
#     return magn, ene


#     def run_ising_dynamics(self):
#         import random as rd

# N = SLRG_obj.system.N
# nmax= 1000
# T = 0.001
# eigV = SLRG_obj.eigV
# H = SLRG_obj.system.H

# def calcEnergy(m, H):
#     '''Energy of a given configuration'''
#     energy = 0
#     for i in range(N):
#         neigh=[w['weight']*m[nn] for nn, w in dict(H[i]).items()]
#         energy += -m[i]*np.sum(neigh)
#     return energy/4.

# #Initialize magnetization
# labels=np.where((eigV[0] < 0))[0]
# lista=['red' if i in labels else 'blue' for i in range(H.number_of_nodes())]
# bigene = []
# for replica in range(1):
#     if replica != 0:
#         m=[1 if i in labels else -1 for i in range(H.number_of_nodes())]
#     else:
#         m = np.random.choice([-1, 1], size=N)

#     #Metropolis
#     magn = []
#     ene = []
#     sample = rd.sample(H.nodes, N)
#     for nsteps in tqdm(range(100)):
#         magn.append(np.sum(m))
#         ene.append(calcEnergy(m, H))
#         for i in range(N):
#             node=sample[i]
#             m_new = -m[node]

#             #Metropolis thing
#             neigh=[w['weight']*m[nn] for nn, w in dict(H[node]).items()]
#             E_old=-m[node]*np.sum(neigh)
#             E_new=-m_new*np.sum(neigh)

#             if E_new<E_old:
#                 m[node]=m_new
#             else:
#                 r=rd.uniform(0, 1)
#                 if (r<np.exp(-(E_new-E_old)/T)):
#                     m[node]=m_new
#     bigene.append(ene)


# def smallest_prob_for_erconn(N, pstart=0.1, halving=0.8, testset=100):
#     p = pstart
#     while True:
#         areallconn = []
#         for _ in range(testset):
#             G = nx.erdos_renyi_graph(N, p, seed=None, directed=False)
#             areallconn.append(nx.is_connected(G))
#         if not all(areallconn):
#             lowest_prob = p
#             break
#         p *= halving
#     return lowest_prob


# def logpsace_prob_erconn(N, pHigh=0.5, pstart=0.1, halving=0.8, testset=100):
#     lowest_prob = smallest_prob_for_erconn(N, pstart, halving, testset)
#     lsp = np.logspace(np.log10(lowest_prob), np.log10(pHigh), num=10)
#     return lsp










# def entropy(G, steps=600, is_signed=False, wTresh=1e-15, t1=-2, t2=5):
#     N = G.number_of_nodes()

#     L, w = get_graph_lspectrum_rw(G, is_signed=is_signed)
#     wSig = w[w > wTresh]

#     t = np.logspace(t1, t2, int(steps))
#     S = np.zeros(len(t))
#     VarL = np.zeros(len(t))

#     for i, tau in enumerate(t):
#         rhoTr = np.exp(-tau * w)
#         Tr = np.nansum(rhoTr)
#         rho = np.divide(rhoTr, Tr)
#         S[i] = -np.nansum(rho * np.log(rho)) / np.log(N)
#         avgrho = np.nansum(np.multiply(w, rhoTr)) / Tr
#         av2rho = np.nansum(np.multiply(np.multiply(w, w), rhoTr)) / Tr
#         VarL[i] = av2rho - avgrho**2

#     dS = np.log(N) * np.diff(1 - S) / np.diff(np.log(t))
#     return 1 - S, dS, VarL, t


# def averaged_Sm1(t1Sm1Avg):
#     lenLs = t1Sm1Avg[0][0].__len__()
#     avgSm1 = np.zeros(lenLs)
#     all_t1 = np.concatenate([i[0] for i in t1Sm1Avg])
#     commonLs = np.logspace(np.log10(min(all_t1)), np.log10(max(all_t1)), lenLs)
#     digitizedLs = []
#     for tS in t1Sm1Avg:
#         digiTmp = np.digitize(tS[0], bins=commonLs) - 1
#         np.add.at(avgSm1, digiTmp, tS[1])
#         digitizedLs.extend(digiTmp)
#     unique, counts = np.unique(digitizedLs, return_counts=True)
#     # if avgSm1.shape != counts.shape:
#     #     counts = np.delete(counts, 0)
#     #     print(avgSm1, unique, counts)
#     try:
#         avgSm1 /= counts
#     except ValueError:
#         print("counts has not the same dimensions of avgSm1.")
#         print(counts.shape, avgSm1.shape)
#     return commonLs, avgSm1


# def radial_correlation(field, neighbor_dict, central_index):
#     # Flatten the field array
#     flattened_field = field.flatten()
#     # Compute the mean of the field values
#     mean_field = np.mean(flattened_field)
#     # Compute the correlation function
#     correlation = []
#     for ord in range(1, int(np.sqrt(len(flattened_field)))):
#         kth_order_neighbors = [n for n, d in neighbor_dict.items() if d == ord]
#         for neighbor_index in kth_order_neighbors:
#             correlations_node = []
#             correlation_node = (
#                 np.mean(
#                     (flattened_field[central_index] - mean_field)
#                     * flattened_field[neighbor_index]
#                 )
#                 / np.std(flattened_field) ** 2
#             )
#             correlations_node.append(correlation_node)
#         correlation.append([ord, np.mean(correlations_node)])

#     return correlation


# # Define a function to calculate the Local Moran's I for a given cell
# def local_moran_i(i, field_array, adj):
#     N = len(field_array)
#     x_avg = np.mean(field_array)
#     x_i = field_array[i]
#     x_j = field_array
#     w_ij = adj[i]
#     m2 = np.sum((x_i - x_avg) ** 2) / N
#     I_i = (x_i - x_avg) / m2 * np.sum(w_ij * (x_j - x_avg))
#     # # Extract the values of the cell and its neighbors
#     # cell_value = field_array[row, col]
#     # cell_neighbors = field_array[max(0, row - wndwS):min(n_rows, row + wndwS + 1),
#     #                  max(0, col - wndwS):min(n_cols, col + wndwS + 1)].flatten()
#     # # Calculate the Moran's I numerator (sum of (x_i - x_mean) * (x_j - x_mean) for all i and j)
#     # numerator = (cell_value - field_mean) * np.sum(cell_neighbors - field_mean)
#     # # Calculate the Moran's I denominator (sum of (x_i - x_mean)^2 for all i)
#     # denominator = np.sum((cell_value - field_mean) ** 2)/(n_rows * n_cols)
#     # # Calculate Moran's I statistic
#     # moran_i = (wndwS + 1) * numerator / denominator
#     return I_i


def ising_spinglass_pmJ_2D_Tcrit(L):
    return L ** (-1.0 / 2)





#
def flip_random_fract_edges(G: Graph, p: float):
    """Flips a fraction p of edges (+1 to -1) of a graph G.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    p : float
       The fraction in [0, 1] of edges to be flipped.
    """
    eset = G.edges()
    nedges = len(eset)
    noflip = int(p * nedges)
    #
    if noflip < 1:
        print("p too small")
        exit()
    rndsmpl = random.sample(range(nedges), noflip)
    #
    # all_weights = {e: 1 for e in eset}
    neg_weights = {e: -1 for i, e in enumerate(eset) if i in rndsmpl}
    #
    nx.set_edge_attributes(G, values=1, name="weight")
    nx.set_edge_attributes(G, values=neg_weights, name="weight")


def eigV_for_lattice2D(side, mode='scipy', howmany=1, **kwargs) -> NDArray:
    l = Lattice2D(side, **kwargs)
    l.flip_random_fract_edges()
    l.compute_k_eigvV(with_routine=mode, k=howmany)
    return l.eigV

def adjust_eigV_for_lattice2D(leigV: NDArray) -> NDArray:
    for i,leV in enumerate(leigV):
        leigV[i] = flip_to_positive_majority(leV)
    return leigV

def eigV_for_lattice2D_ptch(**kwargs) -> NDArray:
    return adjust_eigV_for_lattice2D(eigV_for_lattice2D(**kwargs))

def eigv_for_lattice2D(side, mode: str = "full", **kwargs) -> NDArray:
    l = Lattice2D(side, **kwargs)
    l.flip_random_fract_edges()
    match mode:
        case "full":
            l.compute_laplacian_spectrum_weigV()
        case _ if mode.startswith("some"):
            k = int(mode.split("_")[-1])
            l.compute_k_eigvV(k=k)
    return l.eigv

