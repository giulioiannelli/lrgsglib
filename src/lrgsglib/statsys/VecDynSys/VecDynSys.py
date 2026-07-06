"""Base class for vector-state dynamical systems on signed graphs.

Provides infrastructure for models where each node carries a vector state
(e.g. q-state Potts, XY angle, Heisenberg 3D spin, multi-species).

Hierarchy
---------
::

    DynSys
    +-- VecDynSys (this module)
        +-- PottsModel       (discrete q-state, int32)
        +-- XYModel           (continuous angle, float64)
        +-- HeisenbergModel   (3D unit vector, float64)
        +-- MultiSpeciesModel (k-component discrete, int32)
"""

from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from ..DynSys import DynSys

if TYPE_CHECKING:
    from pathlib import Path

    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph


class VecDynSys(DynSys):
    """Base class for vector-state dynamics on signed graphs.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph substrate.
    q : int
        State dimension (number of states for discrete, vector dimension
        for continuous).
    discrete : bool
        If *True*, states are integers (int32); otherwise float64.
    T : float
        Temperature for Metropolis-based updates.
    **kw
        Forwarded to :class:`DynSys`.
    """

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        q: int = 2,
        discrete: bool = True,
        T: float = 1.0,
        savedisk: bool = True,
        **kw: Any,
    ):
        super().__init__(sg, **kw)
        self.q = q
        self.discrete = discrete
        self.T = T
        # Master switch for on-disk persistence of enabled observables
        # (same semantics as BinDynSys.savedisk): subclasses on the shared
        # observable lifecycle (PottsBase, ...) gate their streams on it;
        # inert for subclasses that do not act on it (legacy PottsModel,
        # XYModel, ...).
        self.savedisk = savedisk
        self._state_dtype = np.int32 if discrete else np.float64
        self._adj: NDArray | None = None

        # Initialise state array with proper shape
        shape = self._state_shape
        self.s: NDArray = np.zeros(shape, dtype=self._state_dtype)

    @property
    def _state_shape(self) -> tuple[int, ...]:
        """Shape depends on model specifics — override in subclass."""
        return (self.N,)

    # ------------------------------------------------------------------
    # Adjacency helper
    # ------------------------------------------------------------------
    def _get_adj_matrix(self) -> NDArray:
        """Return the (cached) signed adjacency matrix as dense float64."""
        if self._adj is None:
            self._adj = self.sg.get_signed_adjacency().astype(np.float64)
        return self._adj

    # ------------------------------------------------------------------
    # State initialisation
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise the vector state.

        Modes
        -----
        uniform / random: Random uniform in [0, q) for discrete,
                          U(0, 2*pi) for continuous.
        zero / zeros:     All zeros.
        custom:           User-supplied array.
        """
        match self.ic:
            case "uniform" | "random" | "rand":
                if self.discrete:
                    self.s = np.random.randint(
                        0, self.q, size=self._state_shape
                    ).astype(np.int32)
                else:
                    self.s = np.random.uniform(
                        0, 2 * np.pi, size=self._state_shape
                    ).astype(np.float64)
            case "zero" | "zeros":
                self.s = np.zeros(self._state_shape, dtype=self._state_dtype)
            case "custom":
                if custom is None:
                    raise ValueError("Provide a custom state array.")
                self.s = np.asarray(custom, dtype=self._state_dtype).copy()
                if self.s.shape != self._state_shape:
                    raise ValueError(
                        f"Shape mismatch: {self.s.shape} != {self._state_shape}"
                    )
            case _:
                raise ValueError(f"Unknown IC: {self.ic!r}")

    # ------------------------------------------------------------------
    # State export
    # ------------------------------------------------------------------
    def export_state(self) -> None:
        """Write state to a binary file for C backend."""
        out_suffix = self.run_id or ''
        fname = self.sg.get_p_fname('s', out_suffix=out_suffix)
        self._ensure_dynpath_exists()
        self.sfout = self.dynpath / fname
        self.s.astype(self._state_dtype).tofile(open(self.sfout, 'wb'))
        self.s_0 = self.s.copy()

    # ------------------------------------------------------------------
    # Dynamics stubs
    # ------------------------------------------------------------------
    def ds1step(self, nd: int) -> None:
        raise NotImplementedError("Subclasses must implement ds1step")

    def run(self, **kw: Any) -> None:
        raise NotImplementedError("Subclasses must implement run")

    def run_py(self, **kw: Any) -> None:
        raise NotImplementedError("Subclasses must implement run_py")

    def run_c(self, **kw: Any) -> None:
        raise NotImplementedError(
            "C backend not configured.  Subclass must also inherit CBackendMixin."
        )
