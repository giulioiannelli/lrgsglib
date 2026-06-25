"""Contact process dynamics for signed graphs.

This module provides two flavours of contact-process dynamics:

* :class:`ContactProcessSIR` implements the infection/recovery process
  parameterised by an infection rate ``mu``. Use ``runlang="C0"`` for the
  C backend or ``runlang="py"`` for pure Python.
* :class:`ContactProcessEI` implements excitation-inhibition dynamics driven
  by ``gamma`` and activation choice. Use ``runlang="C1"`` (final output),
  ``"C1D"`` (density tracking), or ``"C1S"`` (snapshots) for the C backend.

Runlang convention: ``C<digit><letters>``
  - digit 0 = SIR, digit 1 = EI
  - letters: D = density, S = snapshots (bare digit = final state only)
  - Update mode selected via ``--update-mode`` flag (naive, cached, frontier,
    adaptive, gillespie)

Examples
--------
Run the infection/recovery process in Python:

>>> cp = ContactProcessSIR(signed_graph, mu=0.2, runlang="py")
>>> cp.init_contact_dynamics()
>>> cp.run(steps=10, tqdm_on=False)

Use the excitation-inhibition C backend with density tracking:

>>> cp = ContactProcessEI(signed_graph, gamma=1.5, runlang="C1D")
>>> cp.init_contact_dynamics()
>>> cp.run(verbose=False)
"""

from __future__ import annotations

import warnings as _warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Literal, cast

import numpy as np
import tqdm
from numba import njit

from .._c_backend import CBackendMixin
from .._observables import ObservableSet, RowTrajectory, ScalarSeries
from .._solver import SolverBackend
from .._solver_engine import get_solver
from ..BinDynSys import BinDynSys
from ...config.const import BIN
from ...utils.tools.chronometer import time_function_accumulate
from .defaults import (
    CP_DENSITY_FBASE,
    CP_OBS_DENSITY,
    CP_OBS_SNAPSHOTS,
    CP_SNAPSHOT_EVERY,
    CP_SNAPSHOT_MAX_BYTES,
    CP_SNAPSHOTS_FBASE,
    CP_SOLVER_NAME,
)

if TYPE_CHECKING:
    from ...graphs.protocols import SignedGraphProtocol as SignedGraph


# ========================================================================
# Contact Process runlang scheme
# ========================================================================
# New: C<digit><letters> where digit={0=SIR, 1=EI}, letters={D=density, S=snapshots}
# Bare digit = final state only.

_CP_SAVE_LETTERS: dict[str, str] = {
    "D": "density",
    "S": "snapshots",
}

# Deprecated legacy codes → (update_mode, output_mode, new_equivalent)
_CP_EI_DEPRECATED: dict[str, tuple[str, str, str]] = {
    "C1A": ("naive",     "snapshots", "C1S"),
    "C1B": ("naive",     "final",     "C1"),
    "C1C": ("naive",     "density",   "C1D"),
    "C1D": ("adaptive",  "density",   "C1D --update-mode adaptive"),
    "C1E": ("cached",    "density",   "C1D --update-mode cached"),
    "C1F": ("gillespie", "density",   "C1D --update-mode gillespie"),
    "C1G": ("cached",    "snapshots", "C1S --update-mode cached"),
    "CU":  ("naive",     "final",     "C1"),
}


def _parse_cp_ei_runlang(code: str) -> tuple[str, str]:
    """Parse CP EI runlang code → (update_mode, output_mode).

    For new-format codes (C1, C1D, C1S), returns default update_mode="naive".
    For deprecated codes (C1A-G, CU), returns legacy mapping with warning.
    """
    upper = code.strip().upper()

    # Handle CU separately (no digit prefix)
    if upper == "CU":
        _warnings.warn(
            "runlang='CU' is deprecated. Use 'C1' instead.",
            DeprecationWarning, stacklevel=4,
        )
        return "naive", "final"

    if not upper.startswith("C1"):
        raise ValueError(f"Expected C1-based EI code, got '{code}'")

    rest = upper[2:]  # after "C1"

    # New-format: letters are all valid save letters
    if rest == "" or all(c in _CP_SAVE_LETTERS for c in rest):
        if rest:
            # Use first letter for output mode (D or S)
            output_mode = _CP_SAVE_LETTERS[rest[0]]
        else:
            output_mode = "final"
        return "naive", output_mode

    # Check deprecated codes
    full_code = f"C1{rest}"
    if full_code in _CP_EI_DEPRECATED:
        update_mode, output_mode, new_equiv = _CP_EI_DEPRECATED[full_code]
        _warnings.warn(
            f"runlang='{code}' is deprecated. Use '{new_equiv}' instead.",
            DeprecationWarning, stacklevel=4,
        )
        return update_mode, output_mode

    raise ValueError(
        f"Unknown CP runlang code: '{code}'. "
        f"Valid: C0 (SIR), C1[D|S] (EI), or legacy codes."
    )


# ========================================================================
# Numba-optimized activation functions for ContactProcessEI Python backend
# ========================================================================

@njit
def _activation_relu(lambda_val: float) -> float:
    """ReLU activation: P = clip(Lambda, 0, 1).

    Parameters
    ----------
    lambda_val : float
        The weighted sum of neighbor states (gamma_eff * sum(w_ij * s[j]))

    Returns
    -------
    float
        Activation probability in [0, 1]
    """
    # Manual clipping for numba scalar compatibility
    if lambda_val < 0.0:
        return 0.0
    elif lambda_val > 1.0:
        return 1.0
    return lambda_val


@njit
def _activation_tanh(lambda_val: float) -> float:
    """Tanh activation: P = (1 + tanh(Lambda)) / 2.

    Maps tanh output from [-1, 1] to [0, 1] for probability interpretation.

    Parameters
    ----------
    lambda_val : float
        The weighted sum of neighbor states (gamma_eff * sum(w_ij * s[j]))

    Returns
    -------
    float
        Activation probability in [0, 1]
    """
    return (1.0 + np.tanh(lambda_val)) / 2.0


class ContactProcessBase(CBackendMixin, BinDynSys):
    """Shared utilities for contact-process dynamics.

    Parameters
    ----------
    sg : SignedGraph
        Target graph for the simulation.
    state_type : {"binary", "bipolar"}, optional
        State encoding passed through to :class:`BinDynSys`.
    **kwargs : Any
        Additional arguments forwarded to :class:`BinDynSys`.
    """

    dyn_UVclass = "contact_process"

    # CBackendMixin configuration
    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "ContactProcess{}"
    _allowed_c_keys: tuple[str, ...] = ()

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        state_type: Literal["binary", "bipolar"] | None = None,
        savedensity: bool = False,
        snapshot_every: int = CP_SNAPSHOT_EVERY,
        **kwargs: Any,
    ) -> None:
        dynpath = getattr(sg, "path_cntct", None)
        if dynpath is None:
            dynpath = getattr(sg, "path_lrgsg", getattr(sg, "path_data", Path.cwd()))

        super().__init__(
            sg,
            dynpath=dynpath,
            state_type=state_type or "binary",
            **kwargs,
        )
        # Disk-backed observables (Phase 2/4 spine): the active-site density
        # series and the spin-configuration trajectory, behind the reusable
        # ObservableSet. ``savedensity`` enables the density series (off by
        # default to preserve legacy behaviour); ``savedyn`` enables snapshots
        # (``s_t``); ``savedisk`` (from BinDynSys, default True) is the master
        # switch deciding RAM vs disk.
        self.savedensity = bool(savedensity)
        self.snapshot_every = max(1, int(snapshot_every))
        self.observables = ObservableSet([
            ScalarSeries(CP_OBS_DENSITY, dtype=np.float64),
            RowTrajectory(CP_OBS_SNAPSHOTS, ncols=sg.N, dtype=np.int8),
        ])
        self.reset_observables()
        self.sini: np.ndarray | None = None
        self.stderr_path: Path | None = None
        self.stderr_fopen = None
        self.cprogram: list[str | Path] = []
        self.out_id: str = self.out_suffix

    # ------------------------------------------------------------------
    # Disk-backed observables (lazy views over the ObservableSet)
    # ------------------------------------------------------------------
    @property
    def s_t(self):
        """Spin-configuration trajectory: an in-RAM list of per-sweep states
        (``savedisk=False``) or a read-only ``(n_rec, N)`` memmap once streamed
        to disk (``savedisk=True`` + ``savedyn``)."""
        return self.observables[CP_OBS_SNAPSHOTS].data

    @s_t.setter
    def s_t(self, value) -> None:
        self.observables[CP_OBS_SNAPSHOTS].set_rows(value)

    @property
    def density(self) -> list:
        """Active-site density series rho(t) = <s> over MCMC time (populated when
        ``savedensity=True``); lazily disk-backed when persisted."""
        return self.observables[CP_OBS_DENSITY].data

    @density.setter
    def density(self, value) -> None:
        self.observables[CP_OBS_DENSITY].set_values(value)

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset cached observables collected during a run."""

        self.observables.reset()
        self._snap_seen = 0

    # ------------------------------------------------------------------
    # Observable recording / persistence (Python backends)
    # ------------------------------------------------------------------
    def _is_py_runlang(self) -> bool:
        """True for the pure-Python loop (the only backend that streams snapshots
        through the ObservableSet); C backends write their own files."""
        return not self.runlang.upper().startswith("C")

    def _record(self) -> None:
        """Sample the enabled observables for the current sweep (Python loop).

        Records *before* the sweep advances ``self.s`` (matching the legacy
        ``s_t.append(self.s.copy())`` cadence: states at sweep times 0..steps-1).
        """
        if self.savedensity:
            self.observables[CP_OBS_DENSITY].append(self.magnetization())
        if self.savedyn:
            snap = self.observables[CP_OBS_SNAPSHOTS]
            if snap.streaming:
                # Streaming to disk: keep every ``snapshot_every``-th sweep.
                if self._snap_seen % self.snapshot_every == 0:
                    snap.stream_row(self.s)
                self._snap_seen += 1
            else:
                snap.append_row(self.s.copy())

    def _guard_snapshot_size(self) -> None:
        """Refuse before writing if the (subsampled) trajectory would blow the cap."""
        n_rec_est = -(-int(self.steps) // self.snapshot_every)   # ceil(steps / every)
        nbytes = n_rec_est * self.N                              # int8 => 1 byte/elem
        if nbytes > CP_SNAPSHOT_MAX_BYTES:
            raise ValueError(
                f"savedyn trajectory would be ~{nbytes / 2**20:.0f} MiB "
                f"({n_rec_est} snapshots x N={self.N}), exceeding "
                f"CP_SNAPSHOT_MAX_BYTES={CP_SNAPSHOT_MAX_BYTES / 2**20:.0f} MiB. "
                f"Increase snapshot_every (currently {self.snapshot_every}), reduce "
                f"steps, or raise CP_SNAPSHOT_MAX_BYTES in ContactProcess/defaults.py."
            )

    def _begin_outputs(self) -> None:
        """Resolve disk paths and open the snapshot stream before a Python run.

        No-op (observables stay in RAM) unless ``savedisk`` is on; only the
        Python backends stream through the ObservableSet (the C subprocess writes
        its own density/snapshot files).
        """
        snap = self.observables[CP_OBS_SNAPSHOTS]
        snap.reset_stream()
        self._snap_seen = 0
        if not self.savedisk or not self._is_py_runlang():
            return
        suf = self.out_suffix
        if self.savedensity:
            self.observables[CP_OBS_DENSITY].set_path(
                self.dynpath / self.sg.get_p_fname(CP_DENSITY_FBASE, suf, ext=BIN))
        if self.savedyn:
            snap.set_path(
                self.dynpath / self.sg.get_p_fname(CP_SNAPSHOTS_FBASE, suf, ext=BIN))
            self._guard_snapshot_size()
            snap.open_stream()

    def _persist_observables(self) -> None:
        """Finalize disk artifacts and switch attributes to disk-backed views."""
        if not self.savedisk or not self._is_py_runlang():
            return
        snap = self.observables[CP_OBS_SNAPSHOTS]
        # Close a streamed trajectory (no-op if never opened).
        snap.close_stream()
        if self.savedyn and snap.path is not None:
            snap.finalize_batch(subsample=self.snapshot_every)
        if self.savedensity:
            self.observables[CP_OBS_DENSITY].persist()

    def _iter_neighbour_data(self, node: int) -> Iterable[tuple[int, float]]:
        return self.sg.get_neighbors_with_weights(node)

    def init_contact_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the state and export data if required."""
        self._check_c_backend_or_fallback()
        self.reset_observables()
        if custom is not None:
            self.ic = "custom"
        self.init_s(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            exName_arg = exName if exName else self.run_id
            # Export uses build_p_fname which handles underscore via join_non_empty
            self.sg._export_edgel_bin(exName=exName_arg)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        """Initialize dynamics if not already done."""
        # Check sini (set at end of init_contact_dynamics) rather than CbaseName
        # because CbaseName has a class-level default from CBackendMixin
        if not hasattr(self, 'sini') or self.sini is None:
            self.init_contact_dynamics()

    def initialize_run_parameters(self, steps: int | None = None, simref: float | None = None) -> None:
        self._set_time_controls(steps=steps, simref=simref)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def contact_sampling(self, tqdm_on: bool) -> None:
        dsNstep = self.dsNstep()
        nodes = np.arange(self.N)
        iterator = tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
        self._begin_outputs()
        for _ in iterator:
            self._record()
            np.random.shuffle(nodes)
            dsNstep(nodes)
        self._persist_observables()

    def run_py(self, tqdm_on: bool = False) -> None:
        self.contact_sampling(tqdm_on)

    # ------------------------------------------------------------------
    # C backend integration (via CBackendMixin)
    # ------------------------------------------------------------------
    def _c_program_suffix(self) -> str:
        """Extract program suffix for ContactProcess binary name.

        Overridden by subclasses: SIR -> "SIR", EI -> "EI".
        """
        raise NotImplementedError("Subclasses must override _c_program_suffix")

    def _get_cleanup_paths(self) -> list[Path | None]:
        """Return paths to clean up after C run."""
        return [getattr(self, 'sfout', None)]

    @staticmethod
    def _validate_activation(activation: str) -> Literal["tanh", "relu"]:
        normalized = activation.lower()
        if normalized not in {"tanh", "relu"}:
            raise ValueError("activation must be either 'tanh' or 'relu'.")
        return cast(Literal["tanh", "relu"], normalized)

    def _dynamics_out_label(self) -> str:
        gamma = getattr(self, "gamma", None)
        if gamma is not None:
            return f"gamma={float(gamma):.4g}"
        mu = getattr(self, "mu", None)
        if mu is not None:
            return f"mu={float(mu):.12g}"
        return ""

    def _build_out_id(self) -> str:
        """Build identifier used by C backends for output files."""
        from ...utils.basic.strings import join_non_empty
        return join_non_empty("_", self._dynamics_out_label(), self.out_suffix, self.run_id)

    def _build_c_arglist_base(self) -> list[str]:
        """Build base argument list common to all ContactSimulator variants."""
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        self.out_id = self._build_out_id()
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            self._get_datdir_arg(),
            syshape,
            self._c_suffix_arg(self.run_id),
            self.out_id,  # out_id already formatted by _build_out_id()
        ]

    def _build_c_arglist(self) -> list[str]:
        raise NotImplementedError("Subclasses must provide C argument lists.")

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def _resolve_backend(self) -> SolverBackend:
        """Map the runlang code to a solver family: a leading ``C`` selects the C
        subprocess (``C0`` SIR / ``C1*`` EI), anything else the Python loop."""
        return (
            SolverBackend.C if self.runlang.upper().startswith("C")
            else SolverBackend.PY
        )

    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        simref: float | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        """Run contact-process dynamics through the shared solver registry.

        Behaviour-preserving front door for the old runlang if/else: the runlang
        resolves to a backend (``py`` Python loop / ``C`` subprocess), the
        registry hands back its solver, and the C subprocess's transient export
        files are cleaned afterwards when ``clean_export`` is set.
        """
        backend = self._resolve_backend()
        solver = get_solver(CP_SOLVER_NAME, backend)
        solver.supports(self)
        self.check_attribute()
        self.initialize_run_parameters(steps, simref)
        solver.execute(self, tqdm_on=tqdm_on, verbose=verbose)
        if backend is SolverBackend.C and clean_export:
            self.remove_run_c_files()
            self.sg.remove_exported_files()


class ContactProcessSIR(ContactProcessBase):
    """Infection-rate contact process (SIR-style) driven by ``mu``.

    Use this class for the standard infection/recovery dynamics. The Python
    backend (``runlang="py"``) mirrors the logic in :meth:`ds1step`, while the
    ``C0`` runlang targets the ``ContactProcessSIR`` executable.
    """

    _allowed_c_keys = ("C0",)

    def _c_program_suffix(self) -> str:
        return "SIR"

    def __init__(
        self,
        sg: "SignedGraph",
        mu: float = 1.0,
        **kwargs: Any,
    ) -> None:
        super().__init__(sg, **kwargs)
        self.mu = float(mu)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def ds1step(self, node: int) -> None:
        neighbours = list(self._iter_neighbour_data(node))
        if self.s[node]:
            rate = self.mu
            for neighbour, weight in neighbours:
                if weight < 0.0 and self.s[neighbour]:
                    rate += -weight
            prob = 1.0 - np.exp(-rate)
            if prob > 0.0 and np.random.random() < prob:
                self.s[node] = np.int8(0)
        else:
            rate = 0.0
            for neighbour, weight in neighbours:
                if weight > 0.0 and self.s[neighbour]:
                    rate += weight
            prob = 1.0 - np.exp(-rate)
            if prob > 0.0 and np.random.random() < prob:
                self.s[node] = np.int8(1)

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        base_args = self._build_c_arglist_base()
        return [
            base_args[0],  # N
            base_args[1],  # p
            f"{self.mu:.12g}",
            f"{self.steps}",
        ] + base_args[2:]


class ContactProcessEI(ContactProcessBase):
    """Excitation-inhibition contact process targeting the ``ContactProcessEI`` binary.

    Parameters
    ----------
    gamma : float
        Excitation strength (rescaled internally by the average degree to match
        the original ``ContactSimulator1*`` interfaces).
    activation : {"tanh", "relu"}, optional
        Non-linearity used by the C kernels; ignored by other backends.
    num_log_samples : int, optional
        Number of log-spaced sample points for density/snapshot output modes.
    update_mode : {"naive", "cached", "frontier", "adaptive", "gillespie"}, optional
        Update strategy for the C backend (default: "naive").
    output_mode : {"final", "density", "snapshots"}, optional
        Output mode for the C backend (default: "final").

    Notes
    -----
    All C backends (``runlang`` starting with ``C1`` or ``CU``) now target
    the unified ``ContactProcessEI`` binary. Legacy runlang values auto-resolve
    to the appropriate update/output mode combination:

    - ``C1``, ``C1B``: naive + final
    - ``C1A``: naive + snapshots
    - ``C1C``: naive + density
    - ``C1D``: adaptive + density
    - ``C1E``: cached + density
    - ``C1F``: gillespie + density
    - ``C1G``: cached + snapshots

    Python dynamics are also available (``runlang="py"``).
    """

    # Validation handled by _parse_cp_ei_runlang() instead of whitelist
    _allowed_c_keys: tuple[str, ...] = ()

    # Valid update modes for unified ContactSimulator
    _valid_update_modes = frozenset({"naive", "cached", "frontier", "adaptive", "gillespie"})
    # Valid output modes for unified ContactSimulator
    _valid_output_modes = frozenset({"final", "density", "snapshots"})

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        gamma: float,
        activation: Literal["tanh", "relu"] = "tanh",
        num_log_samples: int = 1000,
        update_mode: Literal["naive", "cached", "frontier", "adaptive", "gillespie"] = "naive",
        output_mode: Literal["final", "density", "snapshots"] = "final",
        runlang: str | None = None,
        **kwargs: Any,
    ) -> None:
        if runlang is not None:
            kwargs.setdefault("runlang", runlang)
        else:
            kwargs.setdefault("runlang", "C1")
        super().__init__(sg, **kwargs)
        self.gamma = float(gamma)
        k = sg.Ne / sg.N
        self.gamma_eff = self.gamma / k
        self.activation = self._validate_activation(activation)
        self.num_log_samples = int(num_log_samples)

        # Unified backend parameters
        if update_mode not in self._valid_update_modes:
            raise ValueError(f"update_mode must be one of {self._valid_update_modes}")
        if output_mode not in self._valid_output_modes:
            raise ValueError(f"output_mode must be one of {self._valid_output_modes}")
        self.update_mode = update_mode
        self.output_mode = output_mode

        # Lambda caching data structures (initialized in init_contact_dynamics)
        # Numba-compatible contiguous arrays for Python backend
        self._neigh_indices: np.ndarray | None = None    # Flattened neighbor indices
        self._neigh_weights: np.ndarray | None = None    # Flattened neighbor weights
        self._neigh_offsets: np.ndarray | None = None    # Offsets into flattened arrays
        self._reverse_weights: np.ndarray | None = None  # N x N matrix for reverse edge lookup
        self._lambda_arr: np.ndarray | None = None       # Cached lambda values

        # Pre-select activation function (avoid if-else in hot loop)
        if self.activation == 'relu':
            self._activation_fn = _activation_relu
        elif self.activation == 'tanh':
            self._activation_fn = _activation_tanh
        else:
            raise ValueError(f"Unknown activation: {self.activation}")

    def _init_lambda_cache(self) -> None:
        """Initialize lambda array and neighbor structures for Python backend.

        Creates flattened, contiguous arrays for numba compatibility.
        Uses CSR-like format for efficient storage and access.
        """
        N = self.N

        # First pass: collect neighbor data using inherited method
        neigh_list_ragged = []
        neigh_weights_ragged = []

        for i in range(N):
            neighbors = list(self._iter_neighbour_data(i))
            neigh_indices = [n[0] for n in neighbors]
            neigh_weights = [n[1] for n in neighbors]
            neigh_list_ragged.append(neigh_indices)
            neigh_weights_ragged.append(neigh_weights)

        # Flatten into contiguous arrays for numba
        total_edges = sum(len(nl) for nl in neigh_list_ragged)
        self._neigh_indices = np.empty(total_edges, dtype=np.int32)
        self._neigh_weights = np.empty(total_edges, dtype=np.float64)
        self._neigh_offsets = np.empty(N + 1, dtype=np.int32)

        offset = 0
        for i in range(N):
            self._neigh_offsets[i] = offset
            deg = len(neigh_list_ragged[i])
            if deg > 0:
                self._neigh_indices[offset:offset+deg] = neigh_list_ragged[i]
                self._neigh_weights[offset:offset+deg] = neigh_weights_ragged[i]
            offset += deg
        self._neigh_offsets[N] = offset

        # Build reverse edge weight matrix (N x N, sparse but fast lookup)
        # For large graphs (N > 10000), consider using scipy.sparse
        self._reverse_weights = np.zeros((N, N), dtype=np.float64)
        for i in range(N):
            start = self._neigh_offsets[i]
            end = self._neigh_offsets[i + 1]
            for idx in range(start, end):
                j = self._neigh_indices[idx]
                w_ij = self._neigh_weights[idx]
                self._reverse_weights[j, i] = w_ij  # Store weight from i to j

        # Initialize lambda array from current state
        self._lambda_arr = np.zeros(N, dtype=np.float64)
        for i in range(N):
            start = self._neigh_offsets[i]
            end = self._neigh_offsets[i + 1]
            weighted_sum = 0.0
            for idx in range(start, end):
                j = self._neigh_indices[idx]
                w_ij = self._neigh_weights[idx]
                weighted_sum += w_ij * self.s[j]
            self._lambda_arr[i] = self.gamma_eff * weighted_sum

    def init_contact_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialize contact dynamics with lambda caching for Python backend.

        Parameters
        ----------
        custom : Any, optional
            Custom initial state
        exName : str, optional
            Export name for files
        """
        super().init_contact_dynamics(custom, exName)
        # Initialize lambda cache for Python backend
        if not self.runlang.upper().startswith('C'):
            self._init_lambda_cache()

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def ds1step(self, node: int) -> None:
        """Single-step update using cached lambda (Python backend).

        This method is used by the parent class's contact_sampling() when
        the numba-optimized sweep is not called directly.

        Parameters
        ----------
        node : int
            Node index to update
        """
        if self._lambda_arr is None:
            raise RuntimeError("Lambda cache not initialized. Call init_contact_dynamics() first.")

        # Get probability from cached lambda
        P = self._activation_fn(self._lambda_arr[node])

        # Sample new state
        old_state = self.s[node]
        new_state = np.int8(1 if np.random.random() < P else 0)

        # Update if changed
        if new_state != old_state:
            self.s[node] = new_state
            delta = new_state - old_state

            # Update lambda for all neighbors
            start = self._neigh_offsets[node]
            end = self._neigh_offsets[node + 1]
            for idx in range(start, end):
                j = self._neigh_indices[idx]
                w_ji = self._reverse_weights[j, node]
                self._lambda_arr[j] += self.gamma_eff * w_ji * delta

    @staticmethod
    @njit
    def _sweep_ei_lambda_cache(
        state: np.ndarray,
        lambda_arr: np.ndarray,
        neigh_indices: np.ndarray,
        neigh_weights: np.ndarray,
        neigh_offsets: np.ndarray,
        reverse_weights: np.ndarray,
        gamma_eff: float,
        activation_fn,
        N: int
    ) -> None:
        """Numba-optimized sweep for ContactProcessEI with lambda caching.

        Performs one Monte Carlo sweep (N single-node updates) with random node
        selection. Updates lambda values incrementally when states change.

        Parameters
        ----------
        state : np.ndarray[int8]
            Current state configuration (binary: 0 or 1)
        lambda_arr : np.ndarray[float64]
            Cached lambda values (gamma_eff * sum(w_ji * s[j]))
        neigh_indices : np.ndarray[int32]
            Flattened neighbor indices (CSR format)
        neigh_weights : np.ndarray[float64]
            Flattened neighbor weights parallel to neigh_indices
        neigh_offsets : np.ndarray[int32]
            Offsets into flattened arrays for each node (CSR format)
        reverse_weights : np.ndarray[float64]
            N x N matrix of edge weights for O(1) reverse edge lookup
        gamma_eff : float
            Effective gamma (already rescaled by average degree)
        activation_fn : callable
            Pre-selected activation function (_activation_relu or _activation_tanh)
        N : int
            Number of nodes
        """
        for _ in range(N):
            i = np.random.randint(N)

            # Get activation probability from cached lambda
            P = activation_fn(lambda_arr[i])

            # Sample new state
            old_state = state[i]
            new_state = np.int8(1 if np.random.random() < P else 0)

            # Update if changed
            if new_state != old_state:
                state[i] = new_state
                delta = new_state - old_state

                # Update lambda for all neighbors of i
                start = neigh_offsets[i]
                end = neigh_offsets[i + 1]
                for idx in range(start, end):
                    j = neigh_indices[idx]
                    # Find reverse edge weight (j -> i)
                    w_ji = reverse_weights[j, i]
                    lambda_arr[j] += gamma_eff * w_ji * delta

    def make_sweep_fn(self):
        """Vectorized EI sweep over the cached lambda buffers (overrides the
        scalar ``BinDynSys`` default; see :meth:`DynSys.make_sweep_fn`).

        Requires the Python-backend caches built by :meth:`init_contact_dynamics`
        / :meth:`_init_lambda_cache`.
        """
        if self._lambda_arr is None:
            raise RuntimeError("ContactProcessEI lambda cache is not initialized.")
        if (
            self._neigh_indices is None
            or self._neigh_weights is None
            or self._neigh_offsets is None
        ):
            raise RuntimeError("ContactProcessEI neighbour cache is not initialized.")
        if self._reverse_weights is None:
            raise RuntimeError("ContactProcessEI reverse-weights cache is not initialized.")

        def _sweep() -> None:
            self._sweep_ei_lambda_cache(
                self.s,
                self._lambda_arr,
                self._neigh_indices,
                self._neigh_weights,
                self._neigh_offsets,
                self._reverse_weights,
                self.gamma_eff,
                self._activation_fn,
                self.N,
            )

        return _sweep

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _c_program_suffix(self) -> str:
        return "EI"

    def _build_c_arglist(self) -> list[str]:
        base_args = self._build_c_arglist_base()

        # Parse runlang to get update/output modes (handles new + deprecated)
        update_mode, output_mode = _parse_cp_ei_runlang(self.runlang)
        # Explicit params override parsed defaults
        if self.update_mode != "naive":
            update_mode = self.update_mode
        if self.output_mode != "final":
            output_mode = self.output_mode

        # Unified ContactProcessEI argument format:
        # N p gamma steps datdir syshape run_id out_id activation update_mode output_mode [num_samples]
        args = [
            base_args[0],  # N
            base_args[1],  # p
            f"{self.gamma_eff:.12g}",
            f"{self.steps}",
        ] + base_args[2:] + [
            self.activation,
            update_mode,
            output_mode,
        ]
        # Add num_samples for density/snapshots modes
        if output_mode in ("density", "snapshots"):
            args.append(f"{self.num_log_samples}")
        return args

    def _resolve_backend(self) -> SolverBackend:
        """EI accepts only the unified ``C1*`` / ``CU`` C subprocess or ``py``."""
        u = self.runlang.upper()
        if u.startswith("C1") or u == "CU":
            return SolverBackend.C
        if u == "PY":
            return SolverBackend.PY
        raise ValueError(
            f"ContactProcessEI supports C1*, 'CU', or 'py' backends, "
            f"got '{self.runlang}'"
        )

    def run_py(self, tqdm_on: bool = False) -> None:
        """Python backend: the numba-optimized lambda-cache sweep (overrides the
        base scalar ``contact_sampling`` loop).

        Run setup (``check_attribute`` / ``initialize_run_parameters``) is
        performed by the inherited ``run()`` dispatcher before this is invoked.
        """
        iterator = tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
        self._begin_outputs()
        for _ in iterator:
            self._record()
            # numba-optimized sweep (random node selection + incremental lambda)
            self._sweep_ei_lambda_cache(
                self.s,
                self._lambda_arr,
                self._neigh_indices,
                self._neigh_weights,
                self._neigh_offsets,
                self._reverse_weights,
                self.gamma_eff,
                self._activation_fn,
                self.N,
            )
        self._persist_observables()


def ContactProcess(sg: "SignedGraph", *args: Any, **kwargs: Any):
    """
    Factory for contact-process variants.

    Uses ContactProcessEI when gamma/activation/C1/CU backends are requested,
    otherwise falls back to ContactProcessSIR.
    """
    runlang = kwargs.get("runlang")
    if (
        "gamma" in kwargs
        or "activation" in kwargs
        or "num_log_samples" in kwargs
        or "update_mode" in kwargs
        or "output_mode" in kwargs
        or (isinstance(runlang, str) and (
            runlang.upper().startswith("C1") or runlang.upper() == "CU"
        ))
    ):
        return ContactProcessEI(sg, *args, **kwargs)
    return ContactProcessSIR(sg, *args, **kwargs)


# Register ContactProcess's solver backends (py / C) in the shared solver
# registry. Imported after the classes so that importing ContactProcess wires up
# dispatch; no circular import (``_solvers`` does not import these classes).
from . import _solvers  # noqa: E402,F401
