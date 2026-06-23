"""Base classes and helpers for binary dynamical systems."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from numpy.typing import NDArray

from ..DynSys import DynSys

if TYPE_CHECKING:
    from pathlib import Path

    from ...graphs.protocols import SignedGraphProtocol as SignedGraph

StateType = Literal["bipolar", "binary"]


class BinDynSys(DynSys):
    """Base class for binary dynamical systems on signed graphs.

    Nodes take binary states: bipolar {-1, +1} or binary {0, 1}.

    Parameters
    ----------
    state_type : {'bipolar', 'binary'}, optional
        State encoding: ``'bipolar'`` -> {-1, +1} (default), ``'binary'`` -> {0, 1}.
    savemagn : bool, optional
        If True, the sampling backends record the per-spin magnetization series
        ``magn`` (see :meth:`magnetization`). A shared option for *every* binary
        system -- each engine applies it in its own sampler. Other base
        parameters are inherited from :class:`~lrgsglib.statsys.DynSys.DynSys`.

    Directory Management
    --------------------
    When a BinDynSys subclass (IsingDynamics, VoterModel, ContactProcess) is
    instantiated, it automatically creates its dynamics-specific directory
    (ising/, voter/, or contact/) within the graph's data hierarchy.

    Attributes
    ----------
    magn : list[float]
        Per-spin magnetization time series (populated when
        ``savemagn=True``).
    s_t : list[np.ndarray]
        Time series of state configurations.
    dynpath : Path
        Path to the dynamics-specific data directory.
    run_id : str
        Unique identifier for this simulation run.
    """

    _state_dtype = np.int8

    s_t: list[np.ndarray] = []

    def __init__(
        self,
        sg: "SignedGraph",
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
        dynpath: "Path | None" = None,
        *,
        state_type: StateType = "bipolar",
        savemagn: bool = False,
        savedisk: bool = True,
    ):
        super().__init__(
            sg,
            ic=ic,
            field=field,
            runlang=runlang,
            simref=simref,
            steps=steps,
            savedyn=savedyn,
            seed=seed,
            rndStr=rndStr,
            id_string=id_string,
            out_suffix=out_suffix,
            dynpath=dynpath,
        )
        self._set_state_type(state_type)
        self.s: np.ndarray = np.zeros(self.N, dtype=np.int8)
        # Shared binary observable: per-spin magnetization (see magnetization()).
        self.savemagn = savemagn
        # Master switch for on-disk persistence of enabled observables. When
        # True (default) a subclass that implements it (VoterModel) streams/dumps
        # them under ``dynpath`` and exposes lazy, disk-backed attributes; when
        # False the observables stay purely in RAM. Inert for subclasses that do
        # not act on it (e.g. IsingDynamics keeps its own C ``save`` flags).
        self.savedisk = savedisk
        self.magn: list[float] = []

    # ------------------------------------------------------------------
    # State shape (implements abstract property from DynSys)
    # ------------------------------------------------------------------
    @property
    def _state_shape(self) -> tuple[int, ...]:
        return (self.N,)

    # ------------------------------------------------------------------
    # State-type management
    # ------------------------------------------------------------------
    def _set_state_type(self, state_type: StateType) -> None:
        if state_type not in {"bipolar", "binary"}:
            raise ValueError("state_type must be either 'bipolar' or 'binary'.")
        self._state_type: StateType = state_type

    @property
    def state_type(self) -> StateType:
        return self._state_type

    @property
    def inactive_state(self) -> int:
        return -1 if self.state_type == "bipolar" else 0

    @property
    def active_state(self) -> int:
        return 1

    def magnetization(self, s: NDArray | None = None) -> float:
        """Per-spin magnetization ``M = (1/N) sum_i s_i`` of a binary state.

        The intensive (size-comparable) magnetization defined for *every* binary
        system: for bipolar ``{-1, +1}`` states it lies in ``[-1, 1]``; for
        binary ``{0, 1}`` states it equals the active-site density in ``[0, 1]``.

        Parameters
        ----------
        s : ndarray, optional
            State array to measure; defaults to the current state ``self.s``.
        """
        return float(np.mean(self.s if s is None else s))

    def _coerce_custom_state(self, custom: Any) -> np.ndarray:
        # np.array (not asarray) so the state never aliases the caller's
        # array: in-place sweeps would otherwise corrupt the supplied IC.
        arr = np.array(custom, dtype=np.int8)
        if arr.shape != (self.sg.N,):
            raise ValueError("Custom state must match the number of nodes in the graph.")
        valid_values = {self.inactive_state, self.active_state}
        if not set(np.unique(arr)).issubset(valid_values):
            raise ValueError("Custom state contains values incompatible with the chosen state_type.")
        return arr

    # ------------------------------------------------------------------
    # State manipulation
    # ------------------------------------------------------------------
    def flip_spin(self, nd: int) -> None:
        if self.state_type == "bipolar":
            self.s[nd] = np.int8(-self.s[nd])
        else:
            self.s[nd] = np.int8(
                self.active_state if self.s[nd] == self.inactive_state
                else self.inactive_state
            )

    # ------------------------------------------------------------------
    # State initialisation (implements abstract method from DynSys)
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise the binary state array.

        Delegates to :meth:`init_s` for backwards compatibility.
        """
        self.init_s(custom)

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

    # ------------------------------------------------------------------
    # State export (implements abstract method from DynSys)
    # ------------------------------------------------------------------
    def export_state(self) -> None:
        """Write the current state to disk. Delegates to :meth:`export_s_init`."""
        self.export_s_init()

    def export_s_init(self) -> None:
        out_suffix = self.run_id or ''
        fname = self.sg.get_p_fname('s', out_suffix=out_suffix)
        self.sfout = self.dynpath / fname
        self.s.astype('int8').tofile(open(self.sfout, 'wb'))
        self.s_0 = self.s.copy()

    # ------------------------------------------------------------------
    # Dynamics stubs (subclasses override)
    # ------------------------------------------------------------------
    def ds1step(self, nd: int):
        raise NotImplementedError("Subclasses must implement this method")

    def dsNstep(self):
        return np.vectorize(self.ds1step, excluded="self")

    def make_sweep_fn(self):
        """Standard binary MCMC sweep: a random-order pass of ``ds1step`` over
        all ``N`` nodes (see :meth:`DynSys.make_sweep_fn`).

        The node-order buffer and the vectorized update are built once and reused
        each sweep; reshuffling an array of fixed length consumes the RNG stream
        identically to rebuilding it, so the draw order is reproducible.
        """
        nodes = np.arange(self.N, dtype=np.int64)
        dsNstep = self.dsNstep()

        def _sweep() -> None:
            np.random.shuffle(nodes)
            dsNstep(nodes)

        return _sweep

    def run(self, **kw):
        raise NotImplementedError("Subclasses must implement this method")

    def run_py(self, **kw):
        raise NotImplementedError("Subclasses must implement this method")

    def run_c(self, **kw):
        raise NotImplementedError("Subclasses must implement this method")
