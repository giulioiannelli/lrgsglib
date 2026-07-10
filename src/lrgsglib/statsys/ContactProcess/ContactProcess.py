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

from ...config.const import BIN
from .._c_backend import CBackendMixin
from .._csr import build_graph_csr
from .._observables import ObservableSet, RowTrajectory, ScalarSeries
from .._run_host import RunHostMixin
from .._solver import SolverBackend
from ..BinDynSys import BinDynSys
from ._native import load_cp_native
from .defaults import (
    CP_DENSITY_FBASE,
    CP_EI_ACTIVATION_DEFAULT,
    CP_EI_NLS_DEFAULT,
    CP_EI_UPDATE_MODE_DEFAULT,
    CP_OBS_DENSITY,
    CP_OBS_SNAPSHOTS,
    CP_SIR_BETA_DEFAULT,
    CP_SIR_MU_DEFAULT,
    CP_SNAPSHOT_EVERY,
    CP_SNAPSHOT_MAX_BYTES,
    CP_SNAPSHOTS_FBASE,
    CP_SOLVER_NAME,
)

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph


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
    "C1A": ("naive", "snapshots", "C1S"),
    "C1B": ("naive", "final", "C1"),
    "C1C": ("naive", "density", "C1D"),
    "C1D": ("adaptive", "density", "C1D --update-mode adaptive"),
    "C1E": ("cached", "density", "C1D --update-mode cached"),
    "C1F": ("gillespie", "density", "C1D --update-mode gillespie"),
    "C1G": ("cached", "snapshots", "C1S --update-mode cached"),
    "CU": ("naive", "final", "C1"),
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
            DeprecationWarning,
            stacklevel=4,
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
            DeprecationWarning,
            stacklevel=4,
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


class ContactProcessBase(RunHostMixin, CBackendMixin, BinDynSys):
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
    solver_name = CP_SOLVER_NAME

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
            dynpath = getattr(
                sg, "path_lrgsg", getattr(sg, "path_data", Path.cwd())
            )

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
        self.observables = ObservableSet(
            [
                ScalarSeries(CP_OBS_DENSITY, dtype=np.float64),
                RowTrajectory(CP_OBS_SNAPSHOTS, ncols=sg.N, dtype=np.int8),
            ]
        )
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
        n_rec_est = -(
            -int(self.steps) // self.snapshot_every
        )  # ceil(steps / every)
        nbytes = n_rec_est * self.N  # int8 => 1 byte/elem
        if nbytes > CP_SNAPSHOT_MAX_BYTES:
            raise ValueError(
                f"savedyn trajectory would be ~{nbytes / 2**20:.0f} MiB "
                f"({n_rec_est} snapshots x N={self.N}), exceeding "
                f"CP_SNAPSHOT_MAX_BYTES={CP_SNAPSHOT_MAX_BYTES / 2**20:.0f} MiB. "
                f"Increase snapshot_every (currently {self.snapshot_every}), reduce "
                f"steps, or raise CP_SNAPSHOT_MAX_BYTES in ContactProcess/defaults.py."
            )

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C: one directory per run; leaves publish
    # their token schema via _name_tokens)
    # ------------------------------------------------------------------
    def _lang_token(self) -> str:
        """Backend token for the run dirname.

        Resolved from ``runlang`` when ``run()`` has not stashed
        ``_active_backend`` yet: the manual ``init_contact_dynamics()``
        flow exports the C inputs into the rundir BEFORE ``run()``, so the
        name must not shift between the export and the run.
        """
        backend = getattr(self, "_active_backend", None)
        if backend is None:
            backend = self._resolve_backend()
        return backend.value

    def _cfg_model_block(self) -> dict[str, Any]:
        """Constructor arguments shared by every contact-process flavour;
        the leaves extend with their own physics parameters."""
        return {
            "class": type(self).__name__,
            "state_type": self.state_type,
            "savedensity": self.savedensity,
            "snapshot_every": self.snapshot_every,
            "savedyn": self.savedyn,
            "savedisk": self.savedisk,
            "ic": self.ic,
        }

    def _begin_outputs(self) -> None:
        """Resolve the per-run output directory and open the snapshot stream.

        Layout (Phase C): ONE directory per run under ``dynpath``, named by
        the canonical parameter tokens (``_name_tokens``); the in-process
        backends write the SHORT canonical names (``rho.bin`` /
        ``sout.bin``) plus the ``cfg.json`` sidecar inside it, while the C
        subprocess writes its own sprintf-composed names into the SAME
        rundir (redirected via the composed syshape argument). No-op when
        ``savedisk=False`` or when nothing is persisted (lazy filesystem).
        """
        snap = self.observables[CP_OBS_SNAPSHOTS]
        snap.reset_stream()
        self._snap_seen = 0
        if not self.savedisk:
            return
        if not self._is_py_runlang():
            # C subprocess: the exchange targets the rundir stashed when the
            # arglist was built; only the sidecar is written from Python.
            rundir = (
                self._c_rundir
                if self._c_rundir is not None
                else self._run_output_dir()
            )
            if rundir is not None:
                rundir.mkdir(parents=True, exist_ok=True)
                self._write_cfg_sidecar(rundir)
            return
        if not (self.savedensity or self.savedyn):
            return  # nothing persisted -> no directory (lazy filesystem)
        rundir = self._run_output_dir()
        if rundir is None:  # no naming schema -> legacy flat layout
            suf = self.out_suffix
            if self.savedensity:
                self.observables[CP_OBS_DENSITY].set_path(
                    self.dynpath
                    / self.sg.get_p_fname(CP_DENSITY_FBASE, suf, ext=BIN)
                )
            if self.savedyn:
                snap.set_path(
                    self.dynpath
                    / self.sg.get_p_fname(CP_SNAPSHOTS_FBASE, suf, ext=BIN)
                )
                self._guard_snapshot_size()
                snap.open_stream()
            return
        rundir.mkdir(parents=True, exist_ok=True)
        self._write_cfg_sidecar(rundir)
        if self.savedensity:
            self.observables[CP_OBS_DENSITY].set_path(
                rundir / f"{CP_DENSITY_FBASE}{BIN}"
            )
        if self.savedyn:
            snap.set_path(rundir / f"{CP_SNAPSHOTS_FBASE}{BIN}")
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

    def init_contact_dynamics(
        self, custom: Any = None, exName: str = ""
    ) -> None:
        """Initialise the state and export data if required."""
        self._check_c_backend_or_fallback()
        self.reset_observables()
        if custom is not None:
            self.ic = "custom"
        self.init_s(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            # build_cprogram_command -> _build_c_arglist_base creates the
            # per-run directory (Phase-C naming), so the exports below land
            # inside it.
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            exName_arg = exName if exName else self.run_id
            # Export uses build_p_fname which handles underscore via join_non_empty
            self._c_export_graph_inputs(exName=exName_arg)
        self.sini = self.s.copy()

    def export_s_init(self) -> None:
        """Write the initial state where the C binary will read it.

        With a per-run directory active (``self._c_rundir``), the composed
        syshape makes the binary read ``s_<p><run_id>.bin`` from inside the
        rundir, so the export is redirected there (py-export path == C-read
        path); the flat legacy ``dynpath`` location is used otherwise.
        """
        rundir = self._c_rundir
        if rundir is None:
            super().export_s_init()
            return
        rundir.mkdir(parents=True, exist_ok=True)
        self.sfout = rundir / self.sg.get_p_fname(
            "s", out_suffix=self.run_id or ""
        )
        self.s.astype("int8").tofile(self.sfout)
        self.s_0 = self.s.copy()

    def check_attribute(self) -> None:
        """Initialize dynamics if not already done."""
        # Check sini (set at end of init_contact_dynamics) rather than CbaseName
        # because CbaseName has a class-level default from CBackendMixin
        if not hasattr(self, "sini") or self.sini is None:
            self.init_contact_dynamics()

    def initialize_run_parameters(
        self, steps: int | None = None, simref: float | None = None
    ) -> None:
        self._set_time_controls(steps=steps, simref=simref)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """Bare Python sampling loop (no output lifecycle): the run() front
        door (RunHostMixin) opens/persists the outputs around the solver."""
        dsNstep = self.dsNstep()
        nodes = np.arange(self.N)
        iterator = (
            tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
        )
        for _ in iterator:
            self._record()
            np.random.shuffle(nodes)
            dsNstep(nodes)

    def run_py(self, tqdm_on: bool = False) -> None:
        """Direct Python-backend entry (the registry front door is run())."""
        self._begin_outputs()
        try:
            self._sample_py(tqdm_on=tqdm_on)
        finally:
            self._persist_observables()

    # ------------------------------------------------------------------
    # Native pybind backend hooks (the pb solver polymorphs on these)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        """Setup-time capability gate for ``runlang='pb'`` — the solver calls
        this before any output stream opens. The native kernel
        (``_cp_native.cp_sir_sampling``) implements the SIR
        infection/recovery sweep ONLY; every other contact-process flavour
        is a hard capability error (invariant #3), never a silent fallback.
        """
        raise NotImplementedError(
            "runlang='pb': the native contact-process kernel implements the "
            "SIR infection/recovery dynamics only; "
            f"{type(self).__name__} runs on the python or C backend."
        )

    def _run_pb(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """Execute the native kernel (pb-capable leaves override)."""
        raise NotImplementedError(
            f"runlang='pb' is not wired for {type(self).__name__}; use "
            "runlang='py'."
        )

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
        return [getattr(self, "sfout", None)]

    def remove_run_c_files(self, remove_stderr: bool = True) -> None:
        """Remove the transient C exports; drop an empty non-savedisk rundir."""
        super().remove_run_c_files(remove_stderr)
        self._c_prune_rundir()

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

        return join_non_empty(
            "_", self._dynamics_out_label(), self.out_suffix, self.run_id
        )

    def _build_c_arglist_base(self) -> list[str]:
        """Build base argument list common to all ContactSimulator variants.

        The syshape argument carries the per-run dirname (Phase-C naming;
        see ``CBackendMixin._c_syshape_arg``), so every file the binary
        sprintf-composes — the transient ``s_*`` input and the observable
        outputs — lands inside ``self._c_rundir``.
        """
        datdir = self._get_datdir_arg()
        syshape = self._c_syshape_arg()
        self.out_id = self._build_out_id()
        return [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            datdir,
            syshape,
            self._c_suffix_arg(self.run_id),
            self.out_id,  # out_id already formatted by _build_out_id()
        ]

    def _build_c_arglist(self) -> list[str]:
        raise NotImplementedError("Subclasses must provide C argument lists.")

    # ------------------------------------------------------------------
    # Public interface: run() comes from RunHostMixin (the shared skeleton
    # this class's hand-copied version was hoisted into) — resolve the
    # backend, fetch the registered solver, validate, open outputs, execute,
    # persist, clean the C subprocess's transient exports.
    # ------------------------------------------------------------------


class ContactProcessSIR(ContactProcessBase):
    """Infection/recovery contact process driven by an infection rate ``beta``
    and a recovery rate ``mu``.

    Each susceptible (inactive) site catches the infection at rate
    ``beta * sum(weight)`` over its active positive-edge neighbours; each
    infected (active) site recovers spontaneously at rate ``mu`` (negative edges
    to active neighbours add to the recovery rate — frustration accelerates
    healing). ``beta = 1`` reproduces the historical behaviour where the
    infection rate equalled the summed edge weight. The Python backend
    (``runlang="py"``) mirrors :meth:`ds1step`; ``runlang="pb"`` runs the
    same sweep in-process through the native ``_cp_native`` kernel
    (seed-reproducible, ``beta`` folded into the infection couplings); the
    ``C0`` runlang targets the ``ContactProcessSIR`` executable (mu-only —
    ``beta`` is not wired there).
    """

    _allowed_c_keys = ("C0",)

    def _c_program_suffix(self) -> str:
        return "SIR"

    def __init__(
        self,
        sg: "SignedGraph",
        mu: float = CP_SIR_MU_DEFAULT,
        beta: float = CP_SIR_BETA_DEFAULT,
        **kwargs: Any,
    ) -> None:
        super().__init__(sg, **kwargs)
        self.mu = float(mu)
        self.beta = float(beta)

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C)
    # ------------------------------------------------------------------
    def _name_tokens(self) -> list:
        return [
            ("p", float(self.sg.pflip)),
            ("mu", self.mu),
            ("beta", self.beta),
            ("ns", int(self.steps)),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _cfg_model_block(self) -> dict[str, Any]:
        block = super()._cfg_model_block()
        block.update({"mu": self.mu, "beta": self.beta})
        return block

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
                    rate += self.beta * weight
            prob = 1.0 - np.exp(-rate)
            if prob > 0.0 and np.random.random() < prob:
                self.s[node] = np.int8(1)

    # ------------------------------------------------------------------
    # Native pybind backend (cp_sir_sampling kernel; guards below)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        """What the compiled kernel (``cp_native.cpp``) cannot represent is a
        hard capability error (invariant #3): it implements the binary
        (active = 1 / inactive = 0) SIR sweep of the Python reference loop
        and records the density series only."""
        if self.state_type != "binary":
            raise NotImplementedError(
                "runlang='pb': the native SIR kernel encodes states as "
                "binary {0, 1} (the python loop's activity test is not "
                "representable for bipolar states); "
                f"state_type={self.state_type!r} is python-only."
            )
        if self.savedyn:
            raise NotImplementedError(
                "runlang='pb': the native SIR kernel records the density "
                "series and the final state only, not the per-sweep savedyn "
                "snapshot trajectory; use runlang='py'."
            )

    def _run_pb(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """One in-process native run: same record-then-sweep cadence as
        ``_sample_py`` (density sampled BEFORE each sweep, ``steps``
        records). The kernel receives the SIGNED weights; ``beta`` scales
        the infection couplings only (``ds1step`` applies it to positive
        weights), so it is folded into the weight array here and the kernel
        signature stays mu-only."""
        mod = load_cp_native()
        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        nw = np.where(nw > 0.0, self.beta * nw, nw)
        record = bool(self.savedensity)
        s_out, dens = mod.cp_sir_sampling(
            np.ascontiguousarray(self.s, dtype=np.int8),
            ni,
            nw,
            nptr,
            int(self.N),
            float(self.mu),
            int(self.steps),
            int(self.seed),
            record,
        )
        self.s = np.asarray(s_out, dtype=np.int8)
        if record:
            self.observables[CP_OBS_DENSITY].set_values(
                np.asarray(dens, dtype=np.float64).tolist()
            )

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
    _valid_update_modes = frozenset(
        {"naive", "cached", "frontier", "adaptive", "gillespie"}
    )
    # Valid output modes for unified ContactSimulator
    _valid_output_modes = frozenset({"final", "density", "snapshots"})

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        gamma: float,
        activation: Literal["tanh", "relu"] = CP_EI_ACTIVATION_DEFAULT,
        num_log_samples: int = CP_EI_NLS_DEFAULT,
        update_mode: Literal[
            "naive", "cached", "frontier", "adaptive", "gillespie"
        ] = CP_EI_UPDATE_MODE_DEFAULT,
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
            raise ValueError(
                f"update_mode must be one of {self._valid_update_modes}"
            )
        if output_mode not in self._valid_output_modes:
            raise ValueError(
                f"output_mode must be one of {self._valid_output_modes}"
            )
        self.update_mode = update_mode
        self.output_mode = output_mode

        # Lambda caching data structures (initialized in init_contact_dynamics)
        # Numba-compatible contiguous arrays for Python backend
        self._neigh_indices: np.ndarray | None = (
            None  # Flattened neighbor indices
        )
        self._neigh_weights: np.ndarray | None = (
            None  # Flattened neighbor weights
        )
        self._neigh_offsets: np.ndarray | None = (
            None  # Offsets into flattened arrays
        )
        self._reverse_weights: np.ndarray | None = (
            None  # N x N matrix for reverse edge lookup
        )
        self._lambda_arr: np.ndarray | None = None  # Cached lambda values

        # Pre-select activation function (avoid if-else in hot loop)
        if self.activation == "relu":
            self._activation_fn = _activation_relu
        elif self.activation == "tanh":
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
                self._neigh_indices[offset : offset + deg] = neigh_list_ragged[
                    i
                ]
                self._neigh_weights[offset : offset + deg] = (
                    neigh_weights_ragged[i]
                )
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

    def init_contact_dynamics(
        self, custom: Any = None, exName: str = ""
    ) -> None:
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
        if not self.runlang.upper().startswith("C"):
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
            raise RuntimeError(
                "Lambda cache not initialized. Call init_contact_dynamics() first."
            )

        # Get probability from cached lambda
        P = self._activation_fn(self._lambda_arr[node])

        # Sample new state (active/inactive encoding follows state_type:
        # binary {0, 1} or bipolar {-1, +1}). The lambda update below is written
        # in terms of the signed delta, so it is correct for either encoding.
        old_state = self.s[node]
        new_state = np.int8(
            self.active_state if np.random.random() < P else self.inactive_state
        )

        # Update if changed
        if new_state != old_state:
            self.s[node] = new_state
            delta = int(new_state) - int(old_state)

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
        N: int,
        active_val: int,
        inactive_val: int,
    ) -> None:
        """Numba-optimized sweep for ContactProcessEI with lambda caching.

        Performs one Monte Carlo sweep (N single-node updates) with random node
        selection. Updates lambda values incrementally when states change.

        Parameters
        ----------
        state : np.ndarray[int8]
            Current state configuration (binary {0, 1} or bipolar {-1, +1})
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
        active_val, inactive_val : int
            State encoding: the value written for an active vs inactive node
            (binary 1/0 or bipolar +1/-1). The lambda update uses the signed
            delta, so it stays correct for either encoding.
        """
        for _ in range(N):
            i = np.random.randint(N)

            # Get activation probability from cached lambda
            P = activation_fn(lambda_arr[i])

            # Sample new state
            old_state = state[i]
            new_state = np.int8(
                active_val if np.random.random() < P else inactive_val
            )

            # Update if changed
            if new_state != old_state:
                state[i] = new_state
                delta = np.float64(new_state) - np.float64(old_state)

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
            raise RuntimeError(
                "ContactProcessEI lambda cache is not initialized."
            )
        if (
            self._neigh_indices is None
            or self._neigh_weights is None
            or self._neigh_offsets is None
        ):
            raise RuntimeError(
                "ContactProcessEI neighbour cache is not initialized."
            )
        if self._reverse_weights is None:
            raise RuntimeError(
                "ContactProcessEI reverse-weights cache is not initialized."
            )

        active_val = int(self.active_state)
        inactive_val = int(self.inactive_state)

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
                active_val,
                inactive_val,
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
        args = (
            [
                base_args[0],  # N
                base_args[1],  # p
                f"{self.gamma_eff:.12g}",
                f"{self.steps}",
            ]
            + base_args[2:]
            + [
                self.activation,
                update_mode,
                output_mode,
            ]
        )
        # Add num_samples for density/snapshots modes
        if output_mode in ("density", "snapshots"):
            args.append(f"{self.num_log_samples}")
        return args

    def _c_output_mode(self) -> str:
        """Effective output mode of a C run (runlang code, param override)."""
        _, output_mode = _parse_cp_ei_runlang(self.runlang)
        if self.output_mode != "final":
            output_mode = self.output_mode
        return output_mode

    def _effective_output_mode(self) -> str:
        """Resolved output mode for THIS run: the runlang code wins on the
        C backend (``C1D`` implies density), the constructor argument
        elsewhere."""
        upper = self.runlang.upper()
        if upper.startswith("C1") or upper == "CU":
            return self._c_output_mode()
        return self.output_mode

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C)
    # ------------------------------------------------------------------
    def _name_tokens(self) -> list:
        density_on = self._effective_output_mode() == "density"
        return [
            ("p", float(self.sg.pflip)),
            ("gam", self.gamma),
            ("ns", int(self.steps)),
            (
                "nls",
                self.num_log_samples if density_on else None,
                CP_EI_NLS_DEFAULT,
            ),
            ("act", self.activation, CP_EI_ACTIVATION_DEFAULT),
            ("upd", self.update_mode, CP_EI_UPDATE_MODE_DEFAULT),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _cfg_model_block(self) -> dict[str, Any]:
        block = super()._cfg_model_block()
        block.update(
            {
                "gamma": self.gamma,
                "gamma_eff": self.gamma_eff,
                "activation": self.activation,
                "num_log_samples": self.num_log_samples,
                "update_mode": self.update_mode,
                "output_mode": self._effective_output_mode(),
            }
        )
        return block

    def run_cprogram(self, verbose: bool = False) -> None:
        """Execute the C backend, then reconcile its written outputs.

        Mirrors VoterModel: the ``rho_*`` density series the C binary
        writes (output mode ``'density'``) is read back into the density
        observable, so a C run populates ``self.density`` like a Python
        run — previously the file was written but never loaded.
        """
        super().run_cprogram(verbose)
        if self._c_output_mode() != "density":
            return
        # self.out_id (refreshed by _build_c_arglist_base) carries the full
        # C-side label: gamma=<g>[_<out_suffix>][_<run_id>]. Under Phase-C
        # naming the binary wrote inside the per-run directory (composed
        # syshape), so the read-back joins the exact same rundir + name.
        outdir = self._c_rundir if self._c_rundir is not None else self.dynpath
        rho_path = outdir / self.sg.get_p_fname(
            CP_DENSITY_FBASE, self.out_id, ext=BIN
        )
        if rho_path.exists():
            rho = np.fromfile(rho_path, dtype=np.float64)
            self.observables[CP_OBS_DENSITY].set_values(rho.tolist())

    def _resolve_backend(self) -> SolverBackend:
        """EI accepts only the unified ``C1*`` / ``CU`` C subprocess or
        ``py``; ``pb`` resolves to the PB solver so its capability gate
        (``_pb_check_supported``) raises the hard NotImplementedError."""
        u = self.runlang.upper()
        if u.startswith("C1") or u == "CU":
            return SolverBackend.C
        if u == "PY":
            return SolverBackend.PY
        if u.startswith("PB"):
            return SolverBackend.PB
        raise ValueError(
            f"ContactProcessEI supports C1*, 'CU', 'py', or 'pb' backends, "
            f"got '{self.runlang}'"
        )

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """Python backend: the numba-optimized lambda-cache sweep (overrides
        the base scalar loop). Bare loop — the run() front door (RunHostMixin)
        opens/persists the outputs around the solver.

        Run setup (``check_attribute`` / ``initialize_run_parameters``) is
        performed by the inherited ``run()`` dispatcher before this is invoked.
        """
        iterator = (
            tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
        )
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
                int(self.active_state),
                int(self.inactive_state),
            )


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
        or (
            isinstance(runlang, str)
            and (runlang.upper().startswith("C1") or runlang.upper() == "CU")
        )
    ):
        return ContactProcessEI(sg, *args, **kwargs)
    return ContactProcessSIR(sg, *args, **kwargs)


# Register ContactProcess's solver backends (py / C) in the shared solver
# registry. Imported after the classes so that importing ContactProcess wires up
# dispatch; no circular import (``_solvers`` does not import these classes).
from . import _solvers  # noqa: E402,F401
