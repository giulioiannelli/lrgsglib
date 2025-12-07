"""Base classes and helpers for binary dynamical systems."""

from __future__ import annotations

import os
import random
import time
from pathlib import Path
from typing import Any, Iterable, Literal

import numpy as np
from numpy.typing import NDArray

from ..config.const import *
from ..config.funcs import peq_fstr
from ..nx_patches import SignedGraph
from ..utils.basic.strings import generate_random_id, join_non_empty

ZERO_FIELD = lambda N: np.zeros(N, dtype=np.float64)

StateType = Literal["bipolar", "binary"]


class BinDynSys:
    """Base class for binary dynamical systems on signed graphs.
    
    This class provides the foundation for implementing dynamics where nodes
    can take binary states (either bipolar: {-1, +1} or binary: {0, 1}).
    
    Directory Management
    --------------------
    When a BinDynSys subclass (IsingDynamics, VoterModel, ContactProcess) is
    instantiated, it automatically creates its dynamics-specific directory
    (ising/, voter/, or contact/) within the graph's data hierarchy. This
    ensures that dynamics data is only created when actually needed.
    
    Attributes
    ----------
    s_t : list[np.ndarray]
        Time series of state configurations.
    dynpath : Path
        Path to the dynamics-specific data directory.
    run_id : str
        Unique identifier for this simulation run.
    """
    s_t: list[np.ndarray] = []
    dynpath: Path | str = ""
    run_id: str = ""

    def __init__(
        self,
        sg: SignedGraph,
        ic: str = "uniform",
        field: NDArray = None,
        runlang: str = "py",
        simref: float | None = None,
        steps: int | None = None,
        savedyn: bool = False,
        seed: int = None,
        rndStr: bool = False,
        id_string: str = "",
        out_suffix: str = "",
        dynpath: Path = None,
        *,
        state_type: StateType = "bipolar",
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
        self.field = field if field is not None else ZERO_FIELD(self.N)
        self._set_time_controls(simref=simref, steps=steps, allow_default=True)
        if dynpath is None:
            base_dynpath = getattr(self.sg, "path_data", Path.cwd())
            self.dynpath = Path(base_dynpath)
        else:
            self.dynpath = Path(dynpath)
        # Ensure the dynamics directory exists
        self._ensure_dynpath_exists()
        self._set_state_type(state_type)
        self.s: np.ndarray = np.zeros(self.N, dtype=np.int8)
    
    @property
    def N(self) -> int:
        return self.sg.N

    def _ensure_dynpath_exists(self) -> None:
        """
        Ensure the dynamics directory exists. Creates it if necessary.
        This is called during initialization to set up the directory
        structure for dynamics-specific data.
        """
        self.dynpath.mkdir(parents=True, exist_ok=True)

    def _set_state_type(self, state_type: StateType) -> None:
        if state_type not in {"bipolar", "binary"}:
            raise ValueError("state_type must be either 'bipolar' or 'binary'.")
        self._state_type: StateType = state_type

    # ------------------------------------------------------------------
    # Time-control helpers
    # ------------------------------------------------------------------
    def _set_time_controls(
        self,
        *,
        simref: float | None = None,
        steps: int | None = None,
        allow_default: bool = False,
    ) -> None:
        """
        Configure the simulation horizon in Monte Carlo sweeps.

        Each Monte Carlo step corresponds to ``N`` single-site updates
        (Markov time). Users can specify either ``steps`` (number of
        sweeps) or ``simref`` (size-normalised time, where
        ``steps = simref * N``). When both are provided, ``steps`` takes
        precedence. If neither is provided and ``allow_default`` is True,
        a single sweep is used.
        """
        if steps is None and simref is None:
            if not allow_default and hasattr(self, "steps"):
                return
            steps = 1

        resolved_simref = float(simref) if simref is not None else None
        if resolved_simref is not None and resolved_simref <= 0:
            raise ValueError("simref must be positive.")

        resolved_steps = int(steps) if steps is not None else None
        if resolved_steps is not None and resolved_steps < 1:
            raise ValueError("steps must be a positive integer.")

        if resolved_steps is None and resolved_simref is not None:
            resolved_steps = max(1, int(round(resolved_simref * self.N)))
        if resolved_simref is None and resolved_steps is not None:
            resolved_simref = resolved_steps / float(self.N)

        if resolved_steps is None or resolved_simref is None:
            raise ValueError("Unable to resolve simulation steps; provide simref or steps.")

        self.steps = resolved_steps
        self.simref = resolved_simref

    @property
    def simtime(self) -> int:
        """Compatibility alias mapping to :attr:`steps`."""

        return self.steps

    @simtime.setter
    def simtime(self, value: int) -> None:
        self._set_time_controls(steps=value)

    @property
    def total_single_updates(self) -> int:
        """Total number of single-node updates for the configured run."""

        return self.steps * self.N

    @property
    def state_type(self) -> StateType:
        return self._state_type

    @property
    def inactive_state(self) -> int:
        return -1 if self.state_type == "bipolar" else 0

    @property
    def active_state(self) -> int:
        return 1

    def _coerce_custom_state(self, custom: Any) -> np.ndarray:
        arr = np.asarray(custom, dtype=np.int8)
        if arr.shape != (self.sg.N,):
            raise ValueError("Custom state must match the number of nodes in the graph.")
        valid_values = {self.inactive_state, self.active_state}
        if not set(np.unique(arr)).issubset(valid_values):
            raise ValueError("Custom state contains values incompatible with the chosen state_type.")
        return arr
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
        if self.state_type == "bipolar":
            self.s[nd] = np.int8(-self.s[nd])
        else:
            self.s[nd] = np.int8(self.active_state if self.s[nd] == self.inactive_state else self.inactive_state)
    #
    def init_s(self, custom: Any = None) -> None:
        inactive = self.inactive_state
        active = self.active_state
        match self.ic:
            case "uniform" | "random" | "rand":
                self.s = np.random.choice([inactive, active], size=self.sg.N)
            case _ if self.ic.startswith(("gs", "ground_state")):
                idx = int(self.ic.split("_")[-1])
                eig_state = self.sg.get_eigV_bin_check(which=idx).copy()
                if self.state_type == "binary":
                    self.s = np.where(eig_state > 0, active, inactive)
                else:
                    self.s = eig_state
            case "custom":
                if custom is None:
                    raise ValueError("A custom state must be provided for 'custom' initial condition.")
                self.s = self._coerce_custom_state(custom)
            case "homogeneous" | "homo":
                self.s = np.full(self.sg.N, active, dtype=np.int8)
            case "delta":
                self.s = np.full(self.sg.N, inactive, dtype=np.int8)
                self.s[np.random.randint(self.sg.N)] = active
            case _ if self.ic.startswith("mult_rand") or self.ic.startswith("deltas"):
                try:
                    num = int(self.ic.split("_")[-1])
                except (IndexError, ValueError) as exc:
                    raise ValueError("Invalid 'mult_rand' specification.") from exc
                num = max(1, min(num, self.sg.N))
                self.s = np.full(self.sg.N, inactive, dtype=np.int8)
                indices = np.random.choice(self.sg.N, size=num, replace=False)
                self.s[indices] = active
            case _:
                raise ValueError("Invalid initial condition.")
        self.s = np.asarray(self.s, dtype=np.int8)

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
        out_suffix = self.run_id or ''
        fname = self.sg.get_p_fname('s', out_suffix=out_suffix)
        if not out_suffix:
            # ensure a single trailing underscore when run_id is empty (C expects '_')
            fname = Path(fname.stem + '_' + fname.suffix)
        self.sfout = self.dynpath / fname
        self.s.astype('int8').tofile(open(self.sfout, 'wb'))
        self.s_0 = self.s.copy()
    #
    def export_hfield(self):
        out_suffix = self.run_id or ''
        fname = self.sg.get_p_fname('h', out_suffix=out_suffix)
        if not out_suffix:
            fname = Path(fname.stem + '_' + fname.suffix)
        self.hfout = self.dynpath / fname
        self.field.astype('float64').tofile(open(self.hfout, 'wb'))
