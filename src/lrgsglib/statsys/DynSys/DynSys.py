"""Abstract base class for dynamical systems on signed graphs.

This module defines :class:`DynSys`, the root of the dynamics hierarchy.
All dynamical system classes (binary, continuous, vector) inherit from it.

Hierarchy
---------
::

    DynSys (abstract)
    +-- BinDynSys          (int8 states: {-1,+1} or {0,1})
    +-- ContDynSys          (float64 scalar states, ODE integration)
    +-- VecDynSys           (per-node vector states, MC or ODE)
"""

from __future__ import annotations

import os
import random
import time
from abc import ABC, abstractmethod
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from ...config.const import *
from ...config.funcs import peq_fstr
from ...utils.basic.strings import generate_random_id, join_non_empty

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph
    from .._observables import ObservableSet

ZERO_FIELD = lambda N: np.zeros(N, dtype=np.float64)


class DynSys(ABC):
    """Abstract base class for dynamical systems on signed graphs.

    This class provides the foundation shared by all dynamics families:
    time control, random-seed management, directory management, field
    export, and the generic ``set()`` accessor.

    Subclasses must define:

    * ``_state_dtype`` (class-level): numpy dtype for the state array.
    * ``_state_shape`` (property): shape of the per-instance state array.
    * ``init_state()``: initialise the state array ``self.s``.
    * ``export_state()``: serialise the state to disk.
    * ``ds1step(nd)``: single-node update rule.
    * ``run()``, ``run_py()``, ``run_c()``: simulation drivers.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph on which dynamics are simulated.
    ic : str
        Initial-condition mode (interpretation is subclass-specific).
    field : NDArray or None
        External field vector; defaults to zero.
    runlang : str
        Backend selector (``"py"``, ``"C1b"``, etc.).
    simref : float or None
        Size-normalised time (``steps = simref * N``).
    steps : int or None
        Number of Monte Carlo sweeps / integration steps.
    savedyn : bool
        Whether to save the full state trajectory.
    seed : int or None
        RNG seed (auto-generated if *None*).
    rndStr : bool
        Append a random string to the run identifier.
    id_string : str
        User-supplied identifier prefix.
    out_suffix : str
        Suffix appended to output filenames.
    dynpath : Path or None
        Override for the dynamics-specific output directory.
    """

    # --- Subclass must define ---
    _state_dtype: np.dtype = np.float64  # default; BinDynSys overrides to int8

    s_t: list[np.ndarray] = []
    dynpath: Path | str = ""
    run_id: str = ""

    # Disk-backed observables. A model that persists time series populates this
    # with an ``ObservableSet`` in ``__init__`` (see VoterModel); it stays
    # ``None`` for models that keep observables in plain RAM lists. The reusable
    # accumulate/persist/lazy-load mechanics live in ``statsys._observables`` and
    # the codecs in ``statsys._io``; this slot gives consumers (a batch loader,
    # the live-view GUI) one uniform surface to enumerate and read a run's series.
    observables: "ObservableSet | None" = None

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
        dynpath: Path = None,
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
        # dynpath is lazy: created on first write (state/field exports, the
        # _io observable codecs, _get_datdir_arg for C runs) — constructing a
        # dynamics object must not scatter empty dirs on disk.

    # ------------------------------------------------------------------
    # Abstract interface — subclasses MUST implement
    # ------------------------------------------------------------------

    @property
    @abstractmethod
    def _state_shape(self) -> tuple[int, ...]:
        """Shape of the state array (e.g. ``(N,)`` or ``(N, 3)``)."""
        ...

    @abstractmethod
    def init_state(self, custom: Any = None) -> None:
        """Initialise the state array ``self.s``."""
        ...

    @abstractmethod
    def export_state(self) -> None:
        """Write the current state to disk for C backend consumption."""
        ...

    @abstractmethod
    def ds1step(self, nd: int) -> None:
        """Perform one single-node update at node *nd*."""
        ...

    # ------------------------------------------------------------------
    # Frame production (live-view / animation)
    # ------------------------------------------------------------------
    def make_sweep_fn(self) -> Callable[[], None]:
        """Return a zero-argument callable that advances the model one Monte
        Carlo sweep (``N`` single-node updates) in place.

        Consumed by ``statsys._frames.iter_frames`` to step a model while
        snapshotting its state. Built once per run so per-sweep setup (node
        ordering buffer, vectorized update closure, kernel caches) is hoisted
        out of the loop. Binary models inherit a default from ``BinDynSys``;
        models with a vectorized sweep kernel override it.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not support frame production "
            "(no make_sweep_fn); only implemented for binary dynamics so far."
        )

    @abstractmethod
    def run(self, **kw: Any) -> None:
        """Run the full simulation."""
        ...

    @abstractmethod
    def run_py(self, **kw: Any) -> None:
        """Run the simulation using the Python backend."""
        ...

    @abstractmethod
    def run_c(self, **kw: Any) -> None:
        """Run the simulation using the C backend."""
        ...

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def N(self) -> int:
        return self.sg.N

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
        """Configure the simulation horizon in Monte Carlo sweeps.

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
            raise ValueError(
                "Unable to resolve simulation steps; provide simref or steps."
            )

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

    # ------------------------------------------------------------------
    # Directory management
    # ------------------------------------------------------------------
    def _ensure_dynpath_exists(self) -> None:
        """Create the dynamics output directory if it does not exist."""
        self.dynpath.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Random-seed management
    # ------------------------------------------------------------------
    def __init_randomness__(self, seed: int = None) -> None:
        """Initialise random seeds for ``random`` and ``numpy``.

        If no *seed* is provided, one is generated from the current time
        and the process ID.
        """
        self.seed = seed or (
            (int(time.perf_counter() * 1000) + os.getpid()) % (2**32 - 1)
        )
        random.seed(self.seed)
        np.random.seed(self.seed)
        self.rand_str = generate_random_id()

    # ------------------------------------------------------------------
    # Utility methods
    # ------------------------------------------------------------------
    def set(self, attr_name: str, value: Any) -> None:
        """Set an existing attribute dynamically.

        Raises
        ------
        AttributeError
            If the attribute does not already exist.
        """
        if hasattr(self, attr_name):
            setattr(self, attr_name, value)
        else:
            raise AttributeError(
                f"'{self.__class__.__name__}' "
                f"object has no attribute '{attr_name}'"
            )

    @staticmethod
    def _c_suffix_arg(s: str) -> str:
        """Format suffix for C programs.

        Returns the suffix as-is.  The C code's ``build_str_id()``
        function adds the underscore prefix when constructing filenames.
        """
        return s

    def export_hfield(self) -> None:
        """Write the external field to a binary file."""
        out_suffix = self.run_id or ''
        fname = self.sg.get_p_fname('h', out_suffix=out_suffix)
        self._ensure_dynpath_exists()
        self.hfout = self.dynpath / fname
        self.field.astype('float64').tofile(open(self.hfout, 'wb'))
