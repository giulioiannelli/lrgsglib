import subprocess
import tqdm
#
import numpy as np
#
from typing import Any
#
from .BinDynSys import BinDynSys
#
from ..config.const import *
from ..nx_patches import SignedGraph, get_kth_order_neighbours
from ..utils.basic.strings import generate_random_id, join_non_empty
from ..utils.statsys import boltzmann_factor
from ..utils.tools.chronometer import time_function_accumulate
#
class IsingDynamics(BinDynSys):
    magn = []
    ene = []
    s_t = []
    Ising_clusters = []
    k_B = 1

    def __init__(
        self,
        sg: SignedGraph,
        T: float = 0.,
        *,
        NoClust: int = 1,
        thrmSTEP: int = 20, 
        eqSTEP: int = 20, 
        freq: int = 10,
        nstepsIsing: int = 100,
        save_magnetization: bool = False,
        upd_mode: str = "asynchronous",
        **kwargs
    ) -> None:
        super(IsingDynamics, self).__init__(sg, **kwargs)
        self.T = T
        self.nstepsIsing = nstepsIsing
        self.thrmSTEP = thrmSTEP
        self.eqSTEP = eqSTEP
        self.freq = freq
        self.save_magnetization = save_magnetization
        self.NoClust = NoClust
        self.NeigV = -1
        self.upd_mode = upd_mode
    #
    def neigh_ene(self, neigh: list) -> float:
        return np.sum(neigh) / len(neigh)
    #
    def neigh_wghtmagn(self, node: int,on_g: str = SG_GRAPH_REPR) -> list:
        nd = dict(self.sg.gr[on_g][node])
        return [w["weight"] * self.s[nn] for nn, w in nd.items()]
    #
    def metropolis(self, node):
        neigh = self.neigh_wghtmagn(node)
        neighene = self.neigh_ene(neigh)
        DeltaE = 2 * self.s[node] * neighene
        if DeltaE < 0:
            self.flip_spin(node)
        elif np.random.uniform() < boltzmann_factor(DeltaE, self.T, self.k_B):
            self.flip_spin(node)
    #
    def calc_full_energy(self) -> float:
        neigh_energies = np.array([
            self.neigh_ene(self.neigh_wghtmagn(node)) for node in range(self.sg.N)
        ])
        return -np.dot(self.s, neigh_energies)
    #
    def init_ising_dynamics(self, custom: Any = None, exName: str = ""):
        self.init_s(custom)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            self.export_hfield()
            if self.rndStr and not exName:
                exName = self.run_id            
            self.sg._export_edgel_bin(exName=exName)
        self.sini = self.s.copy()
    #
    def check_attribute(self):
        try:
            getattr(self, f"CbaseName")
        except AttributeError:
            self.init_ising_dynamics()

    def initialize_run_parameters(self, T_ising):
        if T_ising:
            self.T = T_ising

    def build_cprogram_command(self):
        self.CbaseName = f"IsingSimulator{self.runlang[-1]}"
        arglist = [
            f"{self.N}",
            f"{self.T:.3g}",
            f"{self.sg.pflip:.3g}",
            f"{self.NoClust}",
            f"{self.thrmSTEP:.3g}",
            f"{self.eqSTEP}",
            self.sg.path_sgdata.relative_to(Path.cwd()),
            self.sg.syshapePth,
            self.run_id,
            self.out_suffix,
            self.upd_mode,
            f"{self.freq}",
            f"{self.NeigV}"
        ]
        self.cprogram = [LRGSG_CCORE_BIN / self.CbaseName] + arglist
    #
    def setup_stderr_logging(self):
        fname = join_non_empty(
            '_', 
            f"err{self.CbaseName}",
            f"{self.N}",
            f"{self.run_id}",
            f"{self.out_suffix}"
        ) + LOG
        self.stderr_path = LRGSG_LOG / fname
        self.stderr_fopen = open(self.stderr_path, 'w')
    #
    def run_cprogram(self, verbose: bool = False):
        
        result = subprocess.run(self.cprogram, stderr=self.stderr_fopen,
                                stdout=subprocess.PIPE)
        self.s = np.frombuffer(result.stdout, dtype=np.int8)
    #
    def metropolis_sampling(self, tqdm_on):
        metropolis_1step = np.vectorize(self.metropolis, excluded="self")
        if self.save_magnetization:
            def save_magn_array():
                self.s_t.append(self.s)
        else:
            def save_magn_array():
                pass

        sample = list(range(self.sg.N))
        iterator = tqdm(range(self.nstepsIsing)) if tqdm_on \
            else range(self.nstepsIsing)
        self.ene = []
        for _ in iterator:
            self.magn.append(np.sum(self.s))
            self.ene.append(self.calc_full_energy())
            metropolis_1step(sample)
            save_magn_array()
    #
    @time_function_accumulate(auto_log=False)
    def run(self, tqdm_on: bool = True, T_ising: float = None, verbose: bool = False, clean_export: bool = True):
        self.check_attribute()
        self.initialize_run_parameters(T_ising)
        if self.runlang.startswith("C"):
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.metropolis_sampling(tqdm_on)
        

    #
    def find_ising_clusters(self, import_cl: bool = False):
        #can be easily reworked
        if import_cl:
            for i in range(self.NoClust):
                self.Ising_clusters.append(
                    np.fromfile(
                        f"{self.sg.path_ising}cl{i}_{self.sg.std_fname}.bin",
                        dtype=int,
                    )
                )
            self.numIsing_cl = len(self.Ising_clusters)
        if self.Ising_clusters:
            return
        #
        self.sg.compute_k_eigvV(k=self.NoClust)
        eigVbin = self.sg.get_eigV_bin_check_list()
        #
        self.Ising_clusters = []
        for j in range(self.NoClust):
            lnodes = list(self.sg.H.nodes())
            lnodes_tmp = lnodes[:]

            def recursive_search(seed, magn_i, clustertmp):
                neighs = get_kth_order_neighbours(self.sg.H, seed, 1)
                neighs = np.array(
                    [e for e in neighs if e not in set(clustertmp)]
                )
                if not neighs.size:
                    return
                samecluster = np.array(eigVbin[j][neighs] == magn_i)
                if not samecluster.any():
                    return
                neighs_samecluster = list(neighs[samecluster])
                clustertmp.extend(neighs_samecluster)
                for ss in neighs_samecluster:
                    recursive_search(ss, magn_i, clustertmp)

            allclusters = []
            for i in lnodes:
                if i not in lnodes_tmp:
                    continue
                if not lnodes_tmp:
                    break
                #
                clustertmp = []
                clustertmp.extend([i])
                #
                recursive_search(i, eigVbin[j][i], clustertmp)
                lnodes_tmp = [e for e in lnodes_tmp if e not in set(clustertmp)]
                allclusters.append(clustertmp)
            allclusters.sort(key=len, reverse=True)
            self.Ising_clusters.append(allclusters[0])
        self.numIsing_cl = len(self.Ising_clusters)
        if self.runlang.startswith("C"):
            self.sg.export_ising_clust()

    #
    def mapping_nodes_to_clusters(self):
        if not self.Ising_clusters:
            self.find_ising_clusters()
        loc = [x for x in range(len(self.Ising_clusters))]
        self.loc = loc
        node_with_inherclust = [
            [[j, loc[i]] for j in clus]
            for i, clus in enumerate(self.Ising_clusters)
        ]
        self.node_with_inherclust = node_with_inherclust
        node_inherclust_flat = [i for j in node_with_inherclust for i in j]
        self.node_inherclust_flat = node_inherclust_flat
        sorted_list = sorted(node_inherclust_flat, key=lambda x: x[0])
        self.sorted_list = sorted_list
        result_array = np.empty((self.sg.side1, self.sg.side2), dtype=object)
        self.result_array = result_array

        # Fill the result_array with tuples from sorted_list
        for i, sublist in enumerate(sorted_list):
            row, col = divmod(
                i, self.sg.side1
            )  # Calculate the row and column index
            result_array[row, col] = sublist[1]
        self.mapping = result_array

    def _remove_sout(self):
        try:
            self.sfout.unlink()
        except FileNotFoundError:
            pass
    def _remove_hfout(self):
        try:
            self.hfout.unlink()
        except FileNotFoundError:
            pass
    def _remove_stderr(self):
        try:
            self.stderr_path.unlink()
        except FileNotFoundError:
            pass

    def remove_run_c_files(self, remove_stderr: bool = True):
        """
        Remove the output files generated by the C program run.
        
        Parameters
        ----------
        remove_stderr : bool, optional
            If True, removes the stderr file. Default is True.
        """
        self._remove_sout()
        self._remove_hfout()
        if remove_stderr:
            self._remove_stderr()
