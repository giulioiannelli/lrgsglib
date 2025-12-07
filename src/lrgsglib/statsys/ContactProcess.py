"""Contact process dynamics for signed graphs.

This module provides two flavours of contact-process dynamics:

* :class:`ContactProcessSIR` implements the infection/recovery process
  parameterised by an infection rate ``mu``. Use this class with the pure
  Python backend (``runlang="py"``) or the ``C0`` runlang to mirror the
  ``ContactSimulator0`` C kernel.
* :class:`ContactProcessEI` implements excitation-inhibition dynamics driven
  by ``gamma`` and activation choice, mapping to the ``ContactSimulator1*``
  kernels (``runlang`` values ``C1``, ``C1a``, ``C1b``, ``C1c``, ``C1d``,
  ``C1e``). This path encapsulates the degree rescaling expected by the C
  implementations.

Examples
--------
Run the infection/recovery process in Python:

>>> cp = ContactProcessSIR(signed_graph, mu=0.2, runlang="py")
>>> cp.init_contact_dynamics()
>>> cp.run(steps=10, tqdm_on=False)

Use the excitation-inhibition C backend (requires compiled C cores):

>>> cp = ContactProcessEI(signed_graph, gamma=1.5, runlang="C1c")
>>> cp.init_contact_dynamics()
>>> cp.run(verbose=False)
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Iterable, Literal, cast

import numpy as np
import tqdm
from numba import njit

from .BinDynSys import BinDynSys
from ..config.const import LOG, LRGSG_CCORE_BIN, LRGSG_LOG
from ..nx_patches import SignedGraph
from ..utils.basic.strings import join_non_empty
from ..utils.tools.chronometer import time_function_accumulate


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


class ContactProcessBase(BinDynSys):
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
    s_t: list[np.ndarray] = []
    _allowed_c_keys: tuple[str, ...] = ()

    def __init__(
        self,
        sg: SignedGraph,
        *,
        state_type: Literal["binary", "bipolar"] | None = None,
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
        self.reset_observables()
        self.sini: np.ndarray | None = None
        self.stderr_path: Path | None = None
        self.stderr_fopen = None
        self.cprogram: list[str | Path] = []
        self.out_id: str = self.out_suffix

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset cached observables collected during a run."""

        self.s_t = []

    def _iter_neighbour_data(self, node: int) -> Iterable[tuple[int, float]]:
        graph = self.sg.gr[self.sg.on_g]
        if hasattr(graph, "is_directed") and graph.is_directed():
            iterator = graph.predecessors(node)
            for neighbour in iterator:
                data = graph[neighbour][node]
                weight = data.get("weight", 1.0)
                yield neighbour, weight
        else:
            iterator = graph.neighbors(node)
            for neighbour in iterator:
                data = graph[node][neighbour]
                weight = data.get("weight", 1.0)
                yield neighbour, weight

    def init_contact_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the state and export data if required."""

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
            # allow empty exName; SignedGraph._export_edgel_bin will add
            # the required trailing underscore when needed
            self.sg._export_edgel_bin(exName=exName_arg)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        try:
            getattr(self, "CbaseName")
        except AttributeError:
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
        for _ in iterator:
            if self.savedyn:
                self.s_t.append(self.s.copy())
            np.random.shuffle(nodes)
            dsNstep(nodes)

    def run_py(self, tqdm_on: bool = False) -> None:
        self.contact_sampling(tqdm_on)

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _c_program_key(self) -> str:
        if not self.runlang or not self.runlang.upper().startswith("C"):
            raise ValueError("C backend requested but runlang is not a C variant.")
        suffix = self.runlang[1:]
        key = "C0" if suffix == "" else f"C{suffix.upper()}"
        if self._allowed_c_keys and key not in self._allowed_c_keys:
            raise ValueError(f"runlang '{self.runlang}' not supported for {self.__class__.__name__}.")
        return key

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

        return join_non_empty("_", self._dynamics_out_label(), self.out_suffix, self.run_id)

    def _build_c_arglist_base(self) -> list[str]:
        """Build base argument list common to all ContactSimulator variants."""

        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except ValueError:
            datdir = self.sg.path_sgdata
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        self.out_id = self._build_out_id()
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            str(datdir),
            syshape,
            self.run_id,
            self.out_id,
        ]

    def _build_c_arglist(self) -> list[str]:
        raise NotImplementedError("Subclasses must provide C argument lists.")

    def build_cprogram_command(self) -> None:
        c_key = self._c_program_key()
        if c_key == "C0":
            program_suffix = "0"
        elif c_key == "C1":
            program_suffix = "1"
        else:
            program_suffix = c_key[1:].lower()
        self.CbaseName = f"ContactSimulator{program_suffix}"
        arglist = self._build_c_arglist()
        self.cprogram = [LRGSG_CCORE_BIN / self.CbaseName] + arglist

    def setup_stderr_logging(self) -> None:
        fname = join_non_empty(
            '_',
            f"err{self.CbaseName}",
            f"{self.N}",
            self.run_id,
            self.out_suffix,
        ) + LOG
        self.stderr_path = LRGSG_LOG / fname
        self.stderr_path.parent.mkdir(parents=True, exist_ok=True)
        self.stderr_fopen = open(self.stderr_path, 'w')

    def run_cprogram(self, verbose: bool = False) -> None:
        if not self.cprogram:
            raise RuntimeError("C program command has not been initialised.")
        binary_path = Path(self.cprogram[0])
        if not binary_path.exists():
            if self.stderr_fopen and not self.stderr_fopen.closed:
                self.stderr_fopen.close()
            self.stderr_fopen = None
            raise FileNotFoundError(
                f"C backend executable '{binary_path}' was not found. "
                "Build the C components (e.g. via `make c-make`) before running the C backend."
            )
        result = subprocess.run(
            self.cprogram,
            stderr=self.stderr_fopen,
            stdout=subprocess.PIPE,
            check=False,
        )
        if self.stderr_fopen and not self.stderr_fopen.closed:
            self.stderr_fopen.close()
        self.stderr_fopen = None
        state = np.frombuffer(result.stdout, dtype=np.int8)
        if state.size != self.N:
            raise RuntimeError("C backend returned an invalid state configuration.")
        self.s = state.copy()

    def _remove_sfout(self) -> None:
        try:
            self.sfout.unlink()
        except FileNotFoundError:
            pass

    def _remove_stderr(self) -> None:
        if not self.stderr_path:
            return
        try:
            self.stderr_path.unlink()
        except FileNotFoundError:
            pass

    def remove_run_c_files(self, remove_stderr: bool = True) -> None:
        self._remove_sfout()
        if remove_stderr:
            self._remove_stderr()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        simref: float | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        self.check_attribute()
        self.initialize_run_parameters(steps, simref)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.contact_sampling(tqdm_on)


class ContactProcessSIR(ContactProcessBase):
    """Infection-rate contact process (SIR-style) driven by ``mu``.

    Use this class for the standard infection/recovery dynamics. The Python
    backend (``runlang="py"``) mirrors the logic in :meth:`ds1step`, while the
    ``C0`` runlang targets the ``ContactSimulator0`` executable.
    """

    _allowed_c_keys = ("C0",)

    def __init__(
        self,
        sg: SignedGraph,
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
    """Excitation-inhibition contact process targeting ``C1`` kernels.

    Parameters
    ----------
    gamma : float
        Excitation strength (rescaled internally by the average degree to match
        the original ``ContactSimulator1*`` interfaces).
    activation : {"tanh", "relu"}, optional
        Non-linearity used by the C1 kernels; ignored by other backends.
    num_log_samples : int, optional
        Number of log samples used by the ``C1c``, ``C1d``, ``C1e``, ``C1f``, and ``C1g`` variants.

    Notes
    -----
    This class is intended for the C backends (``runlang`` starting with
    ``C1``). Python dynamics are not provided for the excitation-inhibition
    path. The ``C1d`` variant uses adaptive frontier optimization for improved
    performance when density is low (< 0.15). The ``C1e`` variant reuses the
    cached-``lambda`` update scheme without maintaining a frontier and shares
    the ``num_log_samples`` argument with ``C1c`` and ``C1d``. The ``C1f``
    variant introduces a Gillespie-style event-driven loop over the frontier
    for low-density regimes. The ``C1g`` variant combines ``C1e``'s cached-lambda
    updates with configuration snapshots at log-spaced intervals (like ``C1a``).
    """

    _allowed_c_keys = ("C1", "C1A", "C1B", "C1C", "C1D", "C1E", "C1F", "C1G")

    def __init__(
        self,
        sg: SignedGraph,
        *,
        gamma: float,
        activation: Literal["tanh", "relu"] = "tanh",
        num_log_samples: int = 1000,
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

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        base_args = self._build_c_arglist_base()
        key = self._c_program_key()
        args = [
            base_args[0],  # N
            base_args[1],  # p
            f"{self.gamma_eff:.12g}",
            f"{self.steps}",
        ] + base_args[2:] + [
            self.activation,
        ]

        if key == "C1A":
            nSampleLog = getattr(self, "nSampleLog", 100)
            args.append(f"{nSampleLog}")
        elif key in ("C1C", "C1D", "C1E", "C1F", "C1G"):
            args.append(f"{self.num_log_samples}")

        return args

    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        simref: float | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        """Run contact process dynamics (C1* or Python backend).

        Parameters
        ----------
        tqdm_on : bool, optional
            Show progress bar (default: True)
        steps : int, optional
            Number of Monte Carlo sweeps (default: None, uses existing config)
        simref : float, optional
            Size-normalised time (steps = simref * N). Ignored if ``steps`` is provided.
        verbose : bool, optional
            Verbose output (default: False)
        clean_export : bool, optional
            Clean up exported files after run (default: True)
        """
        runlang_upper = self.runlang.upper()

        if runlang_upper.startswith("C1"):
            # Use C backend
            super().run(
                tqdm_on=tqdm_on,
                steps=steps,
                simref=simref,
                verbose=verbose,
                clean_export=clean_export,
            )
        elif runlang_upper == "PY":
            # Use Python backend with numba-optimized sweep
            self.check_attribute()
            self.initialize_run_parameters(steps=steps, simref=simref)

            # Use optimized sweep directly instead of contact_sampling
            iterator = tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
            for _ in iterator:
                if self.savedyn:
                    self.s_t.append(self.s.copy())

                # Call numba-optimized sweep
                self._sweep_ei_lambda_cache(
                    self.s,
                    self._lambda_arr,
                    self._neigh_indices,
                    self._neigh_weights,
                    self._neigh_offsets,
                    self._reverse_weights,
                    self.gamma_eff,
                    self._activation_fn,
                    self.N
                )
        else:
            raise ValueError(
                f"ContactProcessEI supports C1* or 'py' backends, got '{self.runlang}'"
            )


# Backwards compatibility alias
ContactProcess = ContactProcessSIR
