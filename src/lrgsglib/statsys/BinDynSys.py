import random
import time
import os

import numpy as np

from numpy.typing import NDArray
from typing import Any, Iterable

from ..config.const import *
from ..config.funcs import peq_fstr
from ..nx_patches import SignedGraph
from ..utils.basic.strings import generate_random_id, join_non_empty

ZERO_FIELD = lambda N: np.zeros(N, dtype=np.float64)

class BinDynSys:
    s_t = []
    dynpath = ''
    run_id = ''

    def __init__(
        self,
        sg: SignedGraph,
        ic: str = "uniform",
        field: NDArray = None,
        runlang: str = "py",
        simpref: int = 1,
        savedyn: bool = False,
        seed: int = None,
        rndStr: bool = False,
        id_string: str = "",
        out_suffix: str = "",
    ):
        self.__init_randomness__(seed)
        self.sg = sg
        self.ic = ic
        self.runlang = runlang
        self.savedyn = savedyn
        self.rndStr = rndStr
        self.run_id = join_non_empty(
            '_',
            id_string, 
            self.rand_str if rndStr else ''
        )
        self.out_suffix = out_suffix
        self.pflip_id = peq_fstr(self.sg.pflip)
        self.field =  field if field is not None else ZERO_FIELD(self.N)
        self.simtime = simpref * self.N
    
    @property
    def N(self) -> int:
        return self.sg.N
    #
    def __init_randomness__(self, seed: int = None) -> None:
        """
        Initialize the random seed of both the random and numpy libraries.
        If no seed is provided, it will use the current time in milliseconds
        plus the process ID to ensure randomness across different runs.

        Parameters
        ----------
        seed : int, optional
            The seed value to initialize the random number generators. If not
            provided, a seed is generated using the current time and process ID.
        Returns
        -------
        None
        """
        self.seed = seed or ((int(time.perf_counter() * 1000) + os.getpid()) % (2**32 - 1))
        random.seed(self.seed)
        np.random.seed(self.seed)
        #
        self.rand_str = generate_random_id()
        #
    #
    def set(self, attr_name: str, value: Any) -> None:
        """
        Set the value of an attribute dynamically.

        Parameters:
        -----------
        - attr_name (str): The name of the attribute to set.
        - value (Any): The value to assign to the attribute.

        Raises:
        -----------
        - AttributeError: If the attribute does not exist.
        """
        if hasattr(self, attr_name):
            setattr(self, attr_name, value)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' \
                                 object has no attribute '{attr_name}'")
    #
    def flip_spin(self, nd: int) -> None:
        self.s[nd] *= -1
    #
    def init_s(self, custom: Any = None):
        match self.ic:
            case 'uniform'|'random'| 'rand':
                self.s = np.random.choice([-1, 1], size=self.sg.N)
            case _ if self.ic.startswith(("gs", "ground_state")):
                self.n_eigV = int(self.ic.split("_")[-1])
                self.s = self.sg.get_eigV_bin_check(which=self.n_eigV).copy()
            case 'custom':
                self.s = custom
            case 'homogeneous'|'homo':
                self.s = np.ones(self.sg.N)
            case 'delta':
                self.s = -np.ones(self.sg.N)
                self.s[np.random.randint(self.sg.N)] = 1
            case _ if self.ic.startswith('mult_rand'|'deltas'):
                num = int(self.ic.split('_')[-1])
                self.s[np.random.randint(self.sg.N, size=num)]
            case _:
                raise ValueError("Invalid initial condition.")

    def ds1step(self, nd: int):
        raise NotImplementedError("Subclasses must implement this method")

    def dsNstep(self):
        return np.vectorize(self.ds1step, excluded="self")

    def run(self):
        raise NotImplementedError("Subclasses must implement this method")

    def run_py(self):
        raise NotImplementedError("Subclasses must implement this method")

    def run_c(self):
        raise NotImplementedError("Subclasses must implement this method")
    #
    def export_s_init(self):
        self.sfout = self.sg.path_ising / self.sg.get_p_fname('s', self.run_id)
        self.s.astype('int8').tofile(open(self.sfout, 'wb'))
        self.s_0 = self.s.copy()
    #
    def export_hfield(self):
        self.hfout =  self.sg.path_ising /  self.sg.get_p_fname('h', self.run_id)
        self.field.astype('float64').tofile(open(self.hfout, 'wb'))
