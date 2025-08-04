from ..utils.basic.strings import generate_random_id, join_non_empty
from ..utils.tools.chronometer import time_function_accumulate
from .common import *
from .BinDynSys import BinDynSys

class IsingDynamics_DEV(BinDynSys):
    dyn_UVclass = "ising_model"
    id_string_isingdyn = ""

    magn = []
    ene = []
    magnc1 = []
    magn_array_save = []
    Ising_clusters = []
    k_B = 1

    def __init__(
        self,
        sg: SignedGraph = Lattice2D,
        T: float = 1.0,
        runlang: str = "py",
        NoClust: int = 1,
        nstepsIsing: int = 100,
        save_magnetization: bool = False,
        **kwargs
    ) -> None:
        self.T = T
        self.sg = sg
        self.dynpath = self.path_ising
        super(IsingDynamics_DEV, self).__init__(self.sg, **kwargs)
        self.nstepsIsing = nstepsIsing
        self.save_magnetization = save_magnetization
        self.NoClust = NoClust

    def bltzmnn_fact(self, enrgy: float) -> float:
        """
        Calculate the Boltzmann factor for a given energy.

        The Boltzmann factor is given by:
        exp(-enrgy / (k_B * T))

        Parameters:
        -----------
        - enrgy (float): Energy value for which the Boltzmann factor is calculated.

        Returns:
        -----------
        - float: Boltzmann factor.
        """
        return np.exp(-enrgy / (self.k_B * self.T))

    #
    def flip_spin(self, nd: int) -> None:
        """
        Flips the spin of a particle at a specified index within the system.

        Parameters:
        -----------
        - nd (int): The index of the particle/spin to be flipped.

        Returns:
        -----------
        - None

        Modifies the state of the system by flipping the sign of the particle
        at the given index `nd`.
        """
        self.s[nd] *= -1
    #
    def neigh_ene(self, neigh: list) -> float:
        """
        Calculate the average energy of a list of neighboring energy values.

        Parameters:
        -----------
        - neigh (list): A list of neighboring energy values.

        Returns:
        -----------
        - float: The average energy.
        """
        return np.sum(neigh) / len(neigh)
    #
    def neigh_wghtmagn(self, node: int) -> list:
        """
        Calculate the weighted magnitudes of neighboring nodes' values for a given node.

        Parameters:
        -----------
        - node (int): The node for which to calculate the weighted magnitudes of neighbors.

        Returns:
        -----------
        - list: A list of weighted magnitudes.
        """
        node_dict = dict(self.sg.H[node])
        return [w["weight"] * self.s[nn] for nn, w in node_dict.items()]