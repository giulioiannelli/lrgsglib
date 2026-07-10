"""Voter-model dynamics on signed graphs.

Runlang convention (the ``runlang=`` backend selector; short codes mirror
``py`` / ``C<digit>``):
  - ``py`` : pure-Python reference implementation (full rule family).
  - ``pb`` : in-process pybind11 kernel (full rule family, GT-compatible).
  - ``np`` : vectorized NumPy synchronous linear voter.
  - ``cu`` : vectorized CuPy (GPU) synchronous linear voter.
  - ``C<digit><letters>`` : compiled C subprocess; digit 0 = voter dynamics,
    letter ``S`` = file snapshots (a bare ``C0`` writes the final state only).
"""

from __future__ import annotations

import warnings as _warnings
from pathlib import Path
from typing import Any

import numpy as np
import tqdm

from .._c_backend import CBackendMixin
from .._csr import build_graph_csr
from .._run_host import RunHostMixin
from ..BinDynSys import BinDynSys
from ...config.const import BIN, NPZ
from ...utils.statsys import (
    cluster_size_distribution,
    consensus_time_stats as _consensus_time_stats,
    edge_sign_arrays,
    giant_cluster_spans as _giant_cluster_spans,
    interface_density_ensemble as _interface_density_ensemble,
    largest_fraction,
    order_parameter_susceptibility as _order_parameter_susceptibility,
    survival_curve as _survival_curve,
)
from . import _voter_io
from .._observables import (
    HistogramSeries,
    ObservableSet,
    RowTrajectory,
    ScalarSeries,
)
from .defaults import (
    CLUSTER_MODE_CODE,
    CLUSTER_MODES,
    DEFAULT_ABSORBING_EVERY,
    DEFAULT_CLUSTER_MODE,
    DEFAULT_NOISE_EPS,
    DEFAULT_NONLIN_ALPHA,
    DEFAULT_QVOTER_Q,
    DEFAULT_RULE,
    DEFAULT_SOUT_EVERY,
    DEFAULT_UPD_MODE,
    GILLESPIE_RULES,
    LINK_ONLY_RULES,
    PY_ONLY_UPD_MODES,
    RULE_CODE,
    UPD_MODE_CODE,
    VOTER_CLDIST_FBASE,
    VOTER_DEFAULT_OP,
    VOTER_MAGN_FBASE,
    VOTER_OBS_CLDIST,
    VOTER_OBS_MAGN,
    VOTER_OBS_SOUT,
    VOTER_RULES,
    VOTER_SOLVER_NAME,
    VOTER_SOUT_FBASE,
    VOTER_SOUT_MAX_BYTES,
    VOTER_SOUT_NLOG_DEFAULT,
    VOTER_UPD_MODES,
    VOTER_UPD_MODES_PLANNED,
)

# Voter save-mode letters
_VOTER_SAVE_LETTERS: dict[str, str] = {
    "S": "snapshots",
}

# Valid new-style codes (snapshot mode determined by presence of "S")
_VOTER_VALID_CODES: set[str] = {"C0", "C0S"}

# Deprecated legacy codes → (new_equivalent, snapshot_mode)
_VOTER_DEPRECATED: dict[str, tuple[str, bool]] = {
    "C1": ("C0S", True),
}

# Update schedules (Axis B). Phase 3 implements asynchronous, synchronous and
# link updating; gillespie is reserved for Phase 4 and raises
# NotImplementedError so the API never silently lies about the schedule.
_VOTER_MODES_IMPLEMENTED: frozenset[str] = frozenset(VOTER_UPD_MODES)
_VOTER_MODES_PLANNED: frozenset[str] = frozenset(VOTER_UPD_MODES_PLANNED)


class VoterModel(RunHostMixin, CBackendMixin, BinDynSys):
    """Binary voter dynamics with optional C backend.

    In voter dynamics, each node copies the state of a randomly chosen
    neighbor, modified by the edge sign (weight). This creates consensus
    or polarization patterns depending on the graph structure.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph to run dynamics on.
    steps : int, optional
        Number of Monte Carlo sweeps.
    simref : float, optional
        Size-normalized time (steps = simref * N).
    eqSTEP : int, optional
        Legacy alias for steps.
    savemagn : bool, optional
        Inherited from :class:`~lrgsglib.statsys.BinDynSys.BinDynSys`. If True,
        record the per-spin magnetization series ``magn`` (= ``mean(s)``) at
        each sweep; honoured identically across every backend.
    rule : str, optional
        Local update rule (Axis A); default ``'linear'``. One of:

        - ``'linear'`` -- copy ``sign(w_ij) * s_j`` of one uniformly chosen
          neighbour (classic voter).
        - ``'majority'`` -- adopt ``sign(sum_j w_ij s_j)`` with a random
          tie-break (local majority-vote).
        - ``'qvoter'`` -- sample ``q`` neighbours with replacement; copy them
          if their signed opinions are unanimous, else flip with probability
          ``eps`` (non-linear q-voter).
        - ``'nonlinear'`` -- nonlinear (power-alpha) voter: adopt signed
          opinion ``sigma`` with probability ``f_sigma**alpha / (f_+**alpha +
          f_-**alpha)``, where ``f_sigma`` is the fraction of neighbours
          holding ``sigma``. alpha=1 -> linear voter, alpha>1 -> local
          majority, alpha<1 -> minority. Yang et al. 2012, Eq. 1 (see
          ``defaults.py`` ref [5]; they count over agent i + neighbours,
          we use neighbour fractions -- identical power form).
    q : int, optional
        Number of neighbours sampled with replacement for the ``'qvoter'``
        rule (must be >= 1); default 2. Ignored by the other rules.
    eps : float, optional
        Flip probability applied when a ``'qvoter'`` sample is not unanimous
        (must be in ``[0, 1]``); default 0.0. Ignored by the other rules.
    alpha : float, optional
        Nonlinearity exponent for the ``'nonlinear'`` rule (must be >= 0);
        default 1.0. ``alpha=1`` recovers the linear voter, ``alpha>1``
        reinforces the local majority. Ignored by the other rules.
    upd_mode : str, optional
        Update schedule (Axis B); default ``'asynchronous'``. One of
        ``'asynchronous'`` (random-sequential, ``N`` single-node updates per
        sweep), ``'synchronous'`` (double-buffer, all nodes from a frozen
        snapshot), ``'link'`` (edge-update; ``rule='linear'`` only), or
        ``'gillespie'`` (rejection-free continuous-time MC; ``rule='linear'``
        only). All schedules run on every backend (Python, C subprocess,
        pybind).
    absorbing_check : bool or None, optional
        Detect and early-stop at an absorbing (zero-frustration)
        configuration, recording the stopping sweep in ``self.absorbed_at``.
        ``None`` (default) auto-selects from the substrate: ``True`` for an
        unsigned graph (it freezes at consensus), ``False`` for a signed graph
        (frustration lets it fluctuate indefinitely). An explicit bool
        overrides.
    absorbing_every : int, optional
        Test the absorbing condition every this many sweeps when
        ``absorbing_check`` is active (clamped to >= 1); default 1.
    track_clusters : bool, optional
        If True, record at each sampled sweep the connected-component
        (cluster) size distribution of the active-edge subgraph, accumulated
        in ``self.cluster_dist``; default False.
    cluster_mode : str, optional
        Edge-activation predicate for the cluster distribution, used only when
        ``track_clusters=True``; default ``'satisfied'``. ``'satisfied'``
        activates edge ``(i, j)`` when ``sign(w_ij) * s_i * s_j = +1``
        (signed-substrate domains); ``'rawspin'`` when ``s_i == s_j`` (edge
        sign ignored).
    freq : int, optional
        Recording frequency.
    nSampleLog : int, optional
        Number of log-spaced samples (for C backend).
    sout_every : int, optional
        Subsampling stride for the ``savedyn`` spin trajectory: keep every
        ``sout_every``-th recorded sweep (clamped to >= 1); default 1 (full
        trajectory). Bounds the trajectory's on-disk/-RAM footprint.
    sout_nlog : int, optional
        Record this many LOG-SPACED ``savedyn`` snapshots instead of the uniform
        ``sout_every`` stride (the in-process analogue of ``nSampleLog``). When
        ``None`` (default), full/uniform recording is used -- but if that would
        exceed ``VOTER_SOUT_MAX_BYTES`` the run auto-falls back to
        ``VOTER_SOUT_NLOG_DEFAULT`` log-spaced snapshots (with a warning) rather
        than raising.
    sout_force_full : bool, optional
        Record the full ``savedyn`` trajectory even when it exceeds
        ``VOTER_SOUT_MAX_BYTES`` (overrides the log-spaced fallback). Default False.
    savedisk : bool, optional
        Persist *enabled* observables to disk under ``dynpath`` and expose them
        as lazy, disk-backed attributes (default True; inherited from
        ``BinDynSys``). Set False for pure in-RAM behaviour (no files).
    **kwargs
        Additional arguments passed to BinDynSys (e.g. ``savemagn``, ``ic``,
        ``out_suffix``, ``rndStr``).

    Attributes
    ----------
    magn, s_t, cluster_dist
        Lazy, disk-backed observable views (see Notes). ``magn`` is the
        per-spin magnetization series; ``s_t`` the spin-configuration
        trajectory (a read-only ``(n_rec, N)`` memmap when on disk); and
        ``cluster_dist`` the per-sweep cluster-size histograms.

    Notes
    -----
    **Persistence (``savedisk=True``, the default).** Each enabled observable is
    written inside ONE per-run directory under ``dynpath``
    (``data/<lattice>/voter/N=<n>/<run dirname>/``, Phase-C naming; the dirname
    is built from the canonical parameter tokens and a ``cfg.json``
    reproducibility sidecar sits next to the files) when a save flag is set --
    ``savemagn`` -> ``m.bin`` (float64), ``savedyn`` -> ``sout.bin``
    (int8 ``(n_rec, N)``), ``track_clusters`` -> ``cldist.npz`` (the C
    subprocess keeps its own ``m_*``/``sout_*``/``cldist_*`` names inside the
    same run directory) -- and the
    matching attribute becomes lazy: ``magn``/``cluster_dist`` cache the small
    series; ``s_t`` memory-maps the big trajectory (never a full-RAM copy). The
    trajectory streams to disk during the Python loop and is subsampled by
    ``sout_every``; its precomputed size is capped by ``VOTER_SOUT_MAX_BYTES``
    (defaults.py) so extensive runs cannot silently spill terabytes. Writes are
    atomic (temp file + ``os.replace``). Load explicitly with ``load_magn`` /
    ``load_sout(mmap=...)`` / ``load_clusters``, and inspect on-disk bytes with
    ``output_sizes``. The ``cldist`` artifact (per-sweep emergent domain-size
    histograms) is DISTINCT from Ising's ``outcl`` (time-averaged moments of
    fixed eigenvector clusters).

    The rule family (Axis A), update schedules (Axis B), absorbing-state
    detection, and cluster-mode semantics -- with their primary-source
    citations (Castellano 2009 social & non-linear q-voter; Iannelli 2026
    signed Laplacian / frustration) -- are documented in
    :mod:`lrgsglib.statsys.VoterModel.defaults`. The signed-substrate
    generalisation (``sigma_j = sign(w_ij) * s_j``) is a project extension,
    not attributable to those papers.

    Examples
    --------
    >>> voter = VoterModel(lattice, steps=100, runlang='py')
    >>> voter.init_voter_dynamics()
    >>> voter.run(tqdm_on=False)
    """

    dyn_UVclass = "voter_model"
    solver_name = VOTER_SOLVER_NAME

    # CBackendMixin configuration
    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "VoterSimulator{}"
    # Validation handled by suffix method
    _allowed_c_keys: tuple[str, ...] = ()

    # ``magn`` / ``s_t`` / ``cluster_dist`` are lazy, disk-backed properties
    # (see the "Disk-backed observables" section below); their backing fields
    # are initialised in ``reset_observables``.
    _voter_snapshot_mode: bool = False
    # Resolved log-spaced savedyn sample sweeps (None = full/uniform), set per run
    # in _begin_outputs; _sout_sample_set is its O(1)-membership companion (py loop).
    _sout_sample_sweeps: "np.ndarray | None" = None
    _sout_sample_set: "set | None" = None

    def __init__(
        self,
        sg,
        *,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
        rule: str = DEFAULT_RULE,
        q: int = DEFAULT_QVOTER_Q,
        eps: float = DEFAULT_NOISE_EPS,
        alpha: float = DEFAULT_NONLIN_ALPHA,
        upd_mode: str = DEFAULT_UPD_MODE,
        absorbing_check: bool | None = None,
        absorbing_every: int = DEFAULT_ABSORBING_EVERY,
        track_clusters: bool = False,
        cluster_mode: str = DEFAULT_CLUSTER_MODE,
        freq: int = 10,
        nSampleLog: int = 100,
        sout_every: int = DEFAULT_SOUT_EVERY,
        sout_nlog: int | None = None,
        sout_force_full: bool = False,
        **kwargs: Any,
    ) -> None:
        dynpath = getattr(sg, "path_voter", None)
        resolved_steps = steps if steps is not None else eqSTEP
        # Build the observable set BEFORE super().__init__: BinDynSys.__init__
        # assigns ``self.magn = []``, which routes through this class's ``magn``
        # property setter and therefore needs ``self.observables`` to exist. The
        # three voter observables, behind the reusable ObservableSet: a
        # magnetization series (m_*.bin), the spin trajectory (sout_*.bin), and
        # the cluster-size-distribution series (cldist_*.npz). The objects own
        # the codec/lazy-load/streaming mechanics; this class owns the policy
        # (gating flags, sout subsampling, the C/pybind-specific persistence).
        self.observables = ObservableSet(
            [
                ScalarSeries(VOTER_OBS_MAGN, dtype=np.float64),
                RowTrajectory(VOTER_OBS_SOUT, ncols=sg.N, dtype=np.int8),
                HistogramSeries(VOTER_OBS_CLDIST),
            ]
        )
        super().__init__(
            sg, dynpath=dynpath, steps=resolved_steps, simref=simref, **kwargs
        )
        # savemagn is inherited from BinDynSys (shared binary option).
        self.rule = self._validate_rule(rule)
        self.q, self.eps, self.alpha = self._validate_rule_params(q, eps, alpha)
        self.upd_mode = self._validate_upd_mode(upd_mode)
        self._validate_rule_mode(self.rule, self.upd_mode)
        # absorbing_check=None auto-selects from the substrate: an UNSIGNED graph
        # (no negative links) reaches a frozen consensus, so early-stop is on; a
        # SIGNED graph can fluctuate forever (frustration), so it is off. An
        # explicit True/False always overrides.
        self.absorbing_check = (
            (self.sg.count_negative_edges() == 0)
            if absorbing_check is None
            else bool(absorbing_check)
        )
        self.absorbing_every = max(1, int(absorbing_every))
        self.absorbed_at: int | None = None
        self.track_clusters = bool(track_clusters)
        self.cluster_mode = self._validate_cluster_mode(cluster_mode)
        # Active-edge neighbour/sign arrays for the per-record cluster recompute
        # (populated by _setup_cluster_recording when track_clusters is on).
        self._cl_idx: list[np.ndarray] | None = None
        self._cl_b: list[np.ndarray] | None = None
        self.freq = freq
        self.nSampleLog = nSampleLog
        # sout (savedyn) subsampling stride: keep every ``sout_every``-th recorded
        # sweep so the trajectory's on-disk/-RAM footprint stays bounded.
        self.sout_every = max(1, int(sout_every))
        # Log-spaced savedyn sampling: ``sout_nlog`` snapshots (None = uniform
        # ``sout_every``). Auto-engaged with VOTER_SOUT_NLOG_DEFAULT points when the
        # full trajectory would exceed VOTER_SOUT_MAX_BYTES, unless
        # ``sout_force_full`` overrides the cap and records the full trajectory.
        self.sout_nlog = None if sout_nlog is None else max(1, int(sout_nlog))
        self.sout_force_full = bool(sout_force_full)
        self.reset_observables()
        self.sini: np.ndarray | None = None
        self.out_id: str = self.out_suffix
        self.magn_path: Path | None = None
        self._edges_cache: list[tuple[int, int, float]] | None = None

    @property
    def eqSTEP(self) -> int:
        """Compatibility alias for the configured number of sweeps."""
        return self.steps

    @eqSTEP.setter
    def eqSTEP(self, value: int) -> None:
        self._set_time_controls(steps=value)

    @staticmethod
    def _validate_upd_mode(upd_mode: str) -> str:
        """Validate the requested update schedule (see ``upd_mode``)."""
        if upd_mode in _VOTER_MODES_IMPLEMENTED:
            return upd_mode
        if upd_mode in _VOTER_MODES_PLANNED:
            raise NotImplementedError(
                f"upd_mode='{upd_mode}' is planned but not yet implemented; "
                f"use one of {sorted(_VOTER_MODES_IMPLEMENTED)}."
            )
        raise ValueError(
            f"Unknown upd_mode='{upd_mode}'. "
            f"Valid: {sorted(_VOTER_MODES_IMPLEMENTED)}."
        )

    @staticmethod
    def _validate_rule(rule: str) -> str:
        """Validate the requested local update rule (see ``rule``)."""
        if rule in VOTER_RULES:
            return rule
        raise ValueError(f"Unknown rule='{rule}'. Valid: {list(VOTER_RULES)}.")

    @staticmethod
    def _validate_cluster_mode(cluster_mode: str) -> str:
        """Validate the cluster-distribution edge-activation mode."""
        if cluster_mode in CLUSTER_MODES:
            return cluster_mode
        raise ValueError(
            f"Unknown cluster_mode='{cluster_mode}'. Valid: {list(CLUSTER_MODES)}."
        )

    @staticmethod
    def _validate_rule_params(
        q: int,
        eps: float,
        alpha: float,
    ) -> tuple[int, float, float]:
        """Validate the q-voter / nonlinear rule parameters."""
        q, eps, alpha = int(q), float(eps), float(alpha)
        if q < 1:
            raise ValueError(f"q must be >= 1 (q-voter), got {q}.")
        if not 0.0 <= eps <= 1.0:
            raise ValueError(f"eps must be in [0, 1], got {eps}.")
        if alpha < 0.0:
            raise ValueError(f"alpha must be >= 0, got {alpha}.")
        return q, eps, alpha

    @staticmethod
    def _validate_rule_mode(rule: str, upd_mode: str) -> None:
        """Reject (rule, upd_mode) combinations that are ill-defined.

        ``link`` (edge-update) and ``gillespie`` (rejection-free CTMC) are both
        intrinsically copy operations, so each is only meaningful for the linear
        voter rule (ref [1] Sec. III.B.3, p. 601).
        """
        if upd_mode == "link" and rule not in LINK_ONLY_RULES:
            raise ValueError(
                f"upd_mode='link' is defined only for rule in "
                f"{sorted(LINK_ONLY_RULES)} (edge-update is a copy operation); "
                f"got rule='{rule}'."
            )
        if upd_mode == "gillespie" and rule not in GILLESPIE_RULES:
            raise ValueError(
                f"upd_mode='gillespie' is defined only for rule in "
                f"{sorted(GILLESPIE_RULES)} (the rejection-free CTMC is a copy "
                f"operation); got rule='{rule}'."
            )

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset the in-RAM observable buffers and all disk-backing state."""
        # Each Observable clears its RAM buffer/cache, disk path, and (for sout)
        # any open stream + recorded-row count; ``s_t`` then lazily memmaps the
        # file once ``sout`` lives on disk and the RAM cache is dropped.
        self.observables.reset()
        self.absorbed_at = None
        # VoterModel-local subsample counter for the streamed sout (this is run
        # policy, not observable state); see _record / _begin_outputs.
        self._sout_seen: int = 0

    # ------------------------------------------------------------------
    # Disk-backed observables (lazy; see the "Flags save by default" contract)
    # ------------------------------------------------------------------
    @property
    def magn(self) -> list[float]:
        """Per-spin magnetization series (loaded from ``m_*.bin`` on demand)."""
        return self.observables[VOTER_OBS_MAGN].data

    @magn.setter
    def magn(self, value) -> None:
        self.observables[VOTER_OBS_MAGN].set_values(value)

    @property
    def sout_sweep_times(self) -> np.ndarray:
        """Sweep index of each recorded ``s_t`` row.

        ``arange(len(s_t))`` for a full/uniform trajectory; the log-spaced sample
        indices when ``savedyn`` subsampling is active (so ``s_t[i]`` is the
        configuration at sweep ``sout_sweep_times[i]`` -- the x-axis for plots).
        """
        n = len(self.s_t)
        if self._sout_sample_sweeps is not None:
            return np.asarray(self._sout_sample_sweeps[:n], dtype=np.int64)
        return np.arange(n, dtype=np.int64)

    @property
    def cluster_dist(self) -> list[dict[int, int]]:
        """Per-sweep cluster-size histograms (loaded from ``cldist_*.npz``)."""
        return self.observables[VOTER_OBS_CLDIST].data

    @cluster_dist.setter
    def cluster_dist(self, value) -> None:
        self.observables[VOTER_OBS_CLDIST].set_values(value)

    @property
    def s_t(self):
        """Spin-configuration trajectory.

        In-RAM when ``savedisk=False``; otherwise a read-only ``(n_rec, N)``
        memmap of ``sout_*.bin`` (never a full-RAM copy). Iterating yields
        per-sweep configurations either way.
        """
        return self.observables[VOTER_OBS_SOUT].data

    @s_t.setter
    def s_t(self, value) -> None:
        self.observables[VOTER_OBS_SOUT].set_rows(value)

    @property
    def _sout_nrec(self) -> int:
        """Rows actually written to ``sout_*.bin`` (back-compat accessor)."""
        return self.observables[VOTER_OBS_SOUT].nrec

    def init_voter_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the spin configuration and export data if required."""
        self._check_c_backend_or_fallback()
        self.reset_observables()
        self.init_s(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            # build_cprogram_command -> _build_c_arglist creates the per-run
            # directory (Phase-C naming), so the exports below land inside it.
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            # Export uses build_p_fname which handles underscore via join_non_empty
            self._c_export_graph_inputs(exName=self.run_id, adj=True)
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
        # Check sini (set at end of init_voter_dynamics) rather than CbaseName
        # because CbaseName has a class-level default from CBackendMixin
        if self.sini is None:
            self.init_voter_dynamics()

    def initialize_run_parameters(
        self,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
    ) -> None:
        chosen_steps = steps if steps is not None else eqSTEP
        self._set_time_controls(steps=chosen_steps, simref=simref)

    # ------------------------------------------------------------------
    # Python dynamics — reference implementation (Axis A rules × Axis B modes)
    # ------------------------------------------------------------------
    def _apply_rule(self, nd: int, s_read: np.ndarray) -> np.int8:
        """Return the new spin for node ``nd`` under ``self.rule``.

        Neighbour states are read from ``s_read`` (the live array for
        asynchronous updating, a frozen snapshot for synchronous), so one rule
        body serves both schedules. The signed substrate enters through the
        effective neighbour opinion ``sigma_j = sign(w_ij) * s_read[j]`` (a
        negative edge gives an anti-copy). See ``defaults.py`` for citations.
        """
        nbrs = self.sg.get_neighbors_with_weights(nd)
        if not nbrs:
            return np.int8(s_read[nd])
        rule = self.rule
        if rule == "linear":
            # Copy a uniformly chosen neighbour (ref [1] Eq. 5).
            j, w = nbrs[np.random.randint(len(nbrs))]
            return np.int8((-1 if w < 0 else 1) * int(s_read[j]))
        if rule == "majority":
            # Local majority-vote: adopt sign(sum_j w_ij s_j) (ref [1] §III.B.1).
            h = 0.0
            for j, w in nbrs:
                h += w * int(s_read[j])
            if h > 0.0:
                return np.int8(1)
            if h < 0.0:
                return np.int8(-1)
            return np.int8(
                1 if np.random.randint(2) else -1
            )  # random tie-break
        cur = int(s_read[nd])
        if rule == "qvoter":
            # Sample q neighbours WITH replacement; copy if unanimous, else flip
            # with prob eps (ref [2] §II, Eq. 1 — repetition is in the paper).
            deg = len(nbrs)
            first: int | None = None
            unanimous = True
            for k in np.random.randint(0, deg, size=self.q):
                j, w = nbrs[k]
                sigma = (-1 if w < 0 else 1) * int(s_read[j])
                if first is None:
                    first = sigma
                elif sigma != first:
                    unanimous = False
                    break
            if unanimous:
                return np.int8(first)
            return np.int8(-cur if np.random.random() < self.eps else cur)
        if rule == "nonlinear":
            # Reinforcement / nonlinear voter: adopt opinion sigma with
            # probability proportional to (signed fraction holding sigma)^alpha,
            #   P(+1) = f_+^alpha / (f_+^alpha + f_-^alpha),
            # where f_+ is the fraction of neighbours with signed opinion +1.
            # alpha=1 recovers the linear voter; alpha>1 reinforces the local
            # majority; alpha<1 favours the minority. Yang et al. 2012 Eq. 1
            # (defaults.py ref [5]); they count over agent i + neighbours, we
            # use neighbour fractions -- the power form is identical.
            deg = len(nbrs)
            nplus = 0
            for j, w in nbrs:
                if ((-1 if w < 0 else 1) * int(s_read[j])) == 1:
                    nplus += 1
            fp = nplus / deg
            fm = 1.0 - fp
            num = fp**self.alpha
            den = num + fm**self.alpha
            pplus = (num / den) if den > 0.0 else 0.5
            return np.int8(1 if np.random.random() < pplus else -1)
        return np.int8(cur)  # unreachable: rule validated in __init__

    def ds1step(self, nd: int) -> None:
        """Asynchronous in-place update of a single node under ``self.rule``."""
        self.s[nd] = self._apply_rule(nd, self.s)

    # -- absorbing-state detection (signed substrate) -----------------
    def _edges(self) -> list[tuple[int, int, float]]:
        """Cached ``(u, v, weight)`` edge list of the substrate."""
        if self._edges_cache is None:
            self._edges_cache = self.sg.get_edges_with_weights()
        return self._edges_cache

    def count_frustrated_edges(self, s: np.ndarray | None = None) -> int:
        """Number of unsatisfied edges (``w_ij s_i s_j < 0``) in ``s``."""
        s = self.s if s is None else s
        return sum(
            1 for u, v, w in self._edges() if w * int(s[u]) * int(s[v]) < 0
        )

    def is_absorbing(self, s: np.ndarray | None = None) -> bool:
        """True iff ``s`` has zero frustrated edges (the dynamics is frozen).

        A frozen voter configuration exists iff the signed graph is balanced
        (``lambda_min(signed Laplacian) ~ 0``); on a frustrated graph it never
        occurs and the dynamics never freezes (ref [3], Properties (ii)-(iii),
        p. 9). On an unsigned graph this reduces to consensus ``|M| = 1``.
        """
        return self.count_frustrated_edges(s) == 0

    # -- domain / interface observables (single-state + trajectory) ---
    def interface_density(self, s: np.ndarray | None = None) -> float:
        """Density of active (inter-domain) links: frustrated edges / total edges.

        ``rho = (# unsatisfied edges) / (# edges)`` in ``[0, 1]`` -- the voter
        "density of active links" (the boundary density between opinion domains).
        ``rho = 0`` is an absorbing (frozen) configuration; on an unsigned graph
        it is the fraction of edges joining opposite opinions. Defaults to the
        current state ``self.s``.
        """
        n_edges = len(self._edges())
        if n_edges == 0:
            return 0.0
        return self.count_frustrated_edges(s) / n_edges

    def interface_density_series(self) -> list[float]:
        """``rho(t)`` over the recorded trajectory ``s_t`` (needs ``savedyn``)."""
        return [self.interface_density(s) for s in self.s_t]

    def largest_domain_fraction(
        self,
        s: np.ndarray | None = None,
        cluster_mode: str | None = None,
    ) -> float:
        """Giant-domain fraction: largest active-edge domain size / ``N``.

        A *domain* is a connected component of the active-edge subgraph (see
        ``cluster_mode``); returns the largest component size / ``N`` in
        ``[1/N, 1]`` -- a standard coarsening order parameter. Defaults to the
        current state ``self.s`` and the instance ``cluster_mode``. For a whole
        trajectory use :meth:`largest_domain_fraction_series`.
        """
        s = self.s if s is None else s
        mode = (
            self.cluster_mode
            if cluster_mode is None
            else self._validate_cluster_mode(cluster_mode)
        )
        idx, signs = self._gillespie_neighbors()
        b = edge_sign_arrays(signs, mode)
        return largest_fraction([s], idx, b)[0]

    def largest_domain_fraction_series(
        self,
        cluster_mode: str | None = None,
    ) -> list[float]:
        """Giant-domain fraction over the recorded ``s_t`` (needs ``savedyn``)."""
        mode = (
            self.cluster_mode
            if cluster_mode is None
            else self._validate_cluster_mode(cluster_mode)
        )
        idx, signs = self._gillespie_neighbors()
        b = edge_sign_arrays(signs, mode)
        return largest_fraction(self.s_t, idx, b)

    def giant_cluster_spans(
        self,
        s: np.ndarray | None = None,
        cluster_mode: str | None = None,
    ) -> bool:
        """Whether the giant active-edge domain spans (wraps) the lattice torus.

        Reshapes the largest-domain mask of ``s`` to the lattice grid and applies
        the torus-wrap percolation test
        (:func:`lrgsglib.utils.statsys.giant_cluster_spans`). Needs a lattice
        substrate exposing ``syshape`` (square / cubic, von Neumann neighbours);
        for a general graph use :meth:`largest_domain_fraction` instead. Defaults
        to the current state ``self.s`` and the instance ``cluster_mode``.
        """
        shape = getattr(self.sg, "syshape", None)
        if shape is None:
            raise NotImplementedError(
                "giant_cluster_spans needs a lattice substrate with `syshape` "
                "(square/cubic); for a general graph use largest_domain_fraction."
            )
        s = self.s if s is None else s
        mode = (
            self.cluster_mode
            if cluster_mode is None
            else self._validate_cluster_mode(cluster_mode)
        )
        idx, signs = self._gillespie_neighbors()
        b = edge_sign_arrays(signs, mode)
        return _giant_cluster_spans(s, idx, b, shape)

    # -- ensemble observables (static; over a list of finished runs) ---
    @staticmethod
    def consensus_time_stats(runs: "Sequence[VoterModel]") -> dict:
        """Exit / consensus-time statistics over an ensemble of finished runs.

        Reads each run's ``absorbed_at`` (``None`` = never absorbed, right-censored
        at its ``steps``) and delegates to
        :func:`lrgsglib.utils.statsys.consensus_time_stats`; the horizon is the
        max ``steps`` across the ensemble.
        """
        runs = list(runs)
        n_steps = max((r.steps for r in runs), default=0)
        return _consensus_time_stats([r.absorbed_at for r in runs], n_steps)

    @staticmethod
    def survival_curve(
        runs: "Sequence[VoterModel]", n_steps: int | None = None
    ) -> np.ndarray:
        """Ensemble survival curve ``S(t)`` (fraction not yet absorbed by sweep t).

        Delegates to :func:`lrgsglib.utils.statsys.survival_curve`; ``n_steps``
        defaults to the max ``steps`` across ``runs``.
        """
        runs = list(runs)
        if n_steps is None:
            n_steps = max((r.steps for r in runs), default=0)
        return _survival_curve([r.absorbed_at for r in runs], n_steps)

    @staticmethod
    def interface_density_ensemble(runs: "Sequence[VoterModel]"):
        """Ensemble-mean interface-density coarsening curve ``rho(t)`` (+ std).

        Each run needs a recorded trajectory (``savedyn=True``); delegates to
        :func:`lrgsglib.utils.statsys.interface_density_ensemble` over each run's
        :meth:`interface_density_series`.
        """
        return _interface_density_ensemble(
            [r.interface_density_series() for r in runs]
        )

    @staticmethod
    def order_parameter_susceptibility(
        runs: "Sequence[VoterModel]",
        op: str = VOTER_DEFAULT_OP,
        n_sites: int | None = None,
    ):
        """Mean and susceptibility ``chi = n_sites * Var(op)`` of an order
        parameter across an ensemble (disorder / seed realisations).

        ``op`` selects the per-run scalar: ``"largest_domain_fraction"`` (giant
        domain / N, default) or ``"interface_density"``. ``n_sites`` defaults to
        the runs' common ``N``. Delegates to
        :func:`lrgsglib.utils.statsys.order_parameter_susceptibility`.
        """
        runs = list(runs)
        if op == "largest_domain_fraction":
            samples = [r.largest_domain_fraction() for r in runs]
        elif op == "interface_density":
            samples = [r.interface_density() for r in runs]
        else:
            raise ValueError(
                f"unknown op={op!r}; expected 'largest_domain_fraction' or "
                "'interface_density'."
            )
        if n_sites is None and runs:
            n_sites = runs[0].N
        return _order_parameter_susceptibility(samples, n_sites=n_sites)

    # -- per-sweep samplers (Axis B) ----------------------------------
    def _record(self) -> None:
        if self.savemagn:
            self.observables[VOTER_OBS_MAGN].append(self.magnetization())
        if self.savedyn:
            sout = self.observables[VOTER_OBS_SOUT]
            # Log-spaced sampling (self._sout_sample_set) takes precedence over the
            # uniform ``sout_every`` stride; either bounds the trajectory footprint.
            sset = self._sout_sample_set
            take = (
                self._sout_seen in sset
                if sset is not None
                else self._sout_seen % self.sout_every == 0
            )
            if take:
                if sout.streaming:
                    sout.stream_row(self.s)  # off-RAM (Python loop + savedisk)
                else:
                    sout.append_row(self.s.copy())  # in-RAM (savedisk=False)
            self._sout_seen += 1
        if self.track_clusters and self._cl_idx is not None:
            self.observables[VOTER_OBS_CLDIST].append(
                dict(
                    cluster_size_distribution(self.s, self._cl_idx, self._cl_b)
                )
            )

    def _should_stop(self, sweep: int) -> bool:
        if not self.absorbing_check or sweep % self.absorbing_every != 0:
            return False
        if self.is_absorbing():
            self.absorbed_at = sweep
            return True
        return False

    def _sweep_async(self) -> None:
        """N single-node updates, nodes drawn WITH replacement (unified w/ C)."""
        s, N = self.s, self.N
        for _ in range(N):
            nd = int(np.random.randint(N))
            s[nd] = self._apply_rule(nd, s)

    def _sweep_sync(self) -> None:
        """Every node updated simultaneously from a frozen snapshot."""
        snap = self.s.copy()
        new = self.s.copy()
        for nd in range(self.N):
            new[nd] = self._apply_rule(nd, snap)
        self.s = new

    def _sweep_link(self) -> None:
        """N edge-copies per sweep (link/edge-update; rule='linear' only)."""
        edges = self._edges()
        m = len(edges)
        if m == 0:
            return
        for _ in range(self.N):
            u, v, w = edges[int(np.random.randint(m))]
            sign = -1 if w < 0 else 1
            if np.random.randint(2):
                self.s[u] = np.int8(sign * int(self.s[v]))
            else:
                self.s[v] = np.int8(sign * int(self.s[u]))

    # -- cluster-size-distribution recording --------------------------
    def _setup_cluster_recording(self) -> None:
        """Cache the active-edge neighbour/sign arrays used to recompute the
        connected-component size distribution at each recorded sweep.

        The distribution is recomputed from scratch per record (a full
        connected-components pass over the active-edge subgraph, O(N + E)); this
        is correct by construction for any update schedule, so all Python
        schedules support ``track_clusters``. Cleared when not tracking.
        """
        if not self.track_clusters:
            self._cl_idx = self._cl_b = None
            return
        idx, signs = self._gillespie_neighbors()
        self._cl_idx = idx
        self._cl_b = edge_sign_arrays(signs, self.cluster_mode)

    # -- rejection-free CTMC sampler (Axis B: gillespie) --------------
    def _gillespie_neighbors(self) -> tuple[list[np.ndarray], list[np.ndarray]]:
        """Per-node neighbour indices and edge signs (ragged, signed substrate).

        ``signs[i][k] = sign(w_{i,nbr})`` so the effective neighbour opinion is
        ``signs[i][k] * s[idx[i][k]]`` -- a negative edge is an anti-copy, the
        same convention used by ``_apply_rule``.
        """
        idx: list[np.ndarray] = []
        signs: list[np.ndarray] = []
        for nd in range(self.N):
            js: list[int] = []
            sg: list[int] = []
            for j, w in self.sg.get_neighbors_with_weights(nd):
                js.append(j)
                sg.append(-1 if w < 0 else 1)
            idx.append(np.asarray(js, dtype=np.int64))
            signs.append(np.asarray(sg, dtype=np.int8))
        return idx, signs

    def _run_gillespie(self, tqdm_on: bool) -> None:
        """Rejection-free continuous-time linear voter on the signed substrate.

        Each node attempts at rate 1 (so unit continuous time == one sweep ==
        ``N`` attempts); given an attempt it copies a uniform signed neighbour,
        which actually flips it iff the chosen edge is frustrated. Hence node
        ``i`` flips as a Poisson process of rate ``r_i = f_i / deg_i`` where
        ``f_i`` is the number of frustrated edges incident to ``i``. Flips are
        sampled directly (node ``i`` w.p. ``r_i / R``, ``dt ~ Exp(R)``,
        ``R = sum_i r_i``), skipping the null events of the rejection method.

        Magnetization / snapshots are recorded at the integer sweep times
        ``0, 1, ..., steps-1`` (the state is piecewise-constant between flips), so
        the ``magn`` series is directly comparable to the per-sweep samplers. The
        process freezes exactly when ``R = 0`` (no frustrated edge), i.e. at an
        absorbing configuration; with ``absorbing_check`` that stops the run and
        records ``absorbed_at``.
        """
        N = self.N
        s = self.s  # int8, mutated in place so _record() / self.s stay in sync
        idx, signs = self._gillespie_neighbors()
        deg = np.array([a.size for a in idx], dtype=np.int64)
        # f_i = number of frustrated edges incident to i (signs[i]*s[nbr] != s[i]).
        f = np.zeros(N, dtype=np.int64)
        for i in range(N):
            if deg[i]:
                f[i] = int(np.count_nonzero(signs[i] * s[idx[i]] != s[i]))
        rate = np.where(deg > 0, f / np.maximum(deg, 1), 0.0)

        self.absorbed_at = None
        t = 0.0  # continuous time, in sweep units
        recorded = 0  # next integer sweep-time still to record
        steps = self.steps
        while recorded < steps:
            R = float(rate.sum())
            if R <= 0.0:
                # Frozen: no flippable node => absorbing configuration.
                if self.absorbing_check:
                    self._record()
                    self.absorbed_at = recorded
                    recorded += 1
                else:
                    while recorded < steps:  # pad the frozen tail
                        self._record()
                        recorded += 1
                break
            dt = -np.log(np.random.random()) / R
            t_next = t + dt
            # State is constant on [t, t_next): emit any integer times in between.
            while recorded < steps and recorded < t_next:
                self._record()
                recorded += 1
            if recorded >= steps:
                break
            # Pick the flipping node i with probability r_i / R, then flip it.
            i = int(
                np.searchsorted(
                    np.cumsum(rate), np.random.random() * R, side="right"
                )
            )
            if i >= N:
                i = N - 1
            s[i] = np.int8(-int(s[i]))
            # All edges at i toggle frustration; update i and its neighbours.
            f[i] = deg[i] - f[i]
            rate[i] = f[i] / deg[i] if deg[i] else 0.0
            js = idx[i]
            if js.size:
                now_frustrated = signs[i] * s[js] != s[i]
                f[js] += np.where(now_frustrated, 1, -1).astype(np.int64)
                rate[js] = f[js] / deg[js]
            t = t_next

    def voter_sampling(self, tqdm_on: bool) -> None:
        """Run ``self.steps`` sweeps under the configured rule and schedule.

        Records magnetization / snapshots before each sweep. With
        ``absorbing_check=True`` the run stops early once an absorbing
        (zero-frustration) configuration is reached, recording the sweep index
        in ``self.absorbed_at``.
        """
        self.absorbed_at = None
        self._setup_cluster_recording()
        if self.upd_mode == "gillespie":
            self._run_gillespie(tqdm_on)
            return
        sweep_fn = {
            "asynchronous": self._sweep_async,
            "synchronous": self._sweep_sync,
            "link": self._sweep_link,
        }[self.upd_mode]
        iterator = (
            tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
        )
        for t in iterator:
            self._record()
            if self._should_stop(t):
                break
            sweep_fn()

    def run_py(self, tqdm_on: bool = False) -> None:
        self.voter_sampling(tqdm_on)

    # ------------------------------------------------------------------
    # Native-backend capability guard
    # ------------------------------------------------------------------
    def _assert_native_supports_config(self, backend: str) -> None:
        """Native backends implement the full rule family + sampler axis
        (including the ``gillespie`` CTMC, via the shared ``_ccore`` kernel) +
        absorbing early-stop. The **pybind** backend also records the per-sweep
        ``savedyn`` trajectory in memory (returned as ``s_t``); the **C
        subprocess** does not (use ``C0S`` file snapshots there). Refuse rather
        than silently return an empty ``s_t``. ``PY_ONLY_UPD_MODES`` (currently
        empty) lists any schedule that is still Python-only.
        """
        if self.savedyn and not self.runlang.lower().startswith("pb"):
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not record the "
                f"per-sweep savedyn trajectory; use runlang='pb' (in-memory "
                f"s_t), runlang='py', or 'C0S' (nSampleLog file snapshots)."
            )
        if self.upd_mode in PY_ONLY_UPD_MODES:
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not implement "
                f"upd_mode='{self.upd_mode}' yet (native CTMC kernel pending); "
                f"use runlang='py'."
            )
        # track_clusters is supported by BOTH native backends: pybind builds the
        # cluster_dist in memory; the C subprocess writes the sparse cldist file
        # (read back by run_cprogram via _voter_io.load_cldist_bin).

    # ------------------------------------------------------------------
    # C backend integration (via CBackendMixin)
    # ------------------------------------------------------------------
    def _c_program_suffix(self) -> str:
        """Map runlang to unified VoterSimulator (always returns '').

        Sets ``_voter_snapshot_mode`` as a side effect.

        New codes: C0 (final), C0S (snapshots).
        Deprecated: C1 → C0S.
        """
        upper = self.runlang.strip().upper()
        if upper in _VOTER_VALID_CODES:
            self._voter_snapshot_mode = "S" in upper[2:]
            return ""
        if upper in _VOTER_DEPRECATED:
            new_equiv, snapshot = _VOTER_DEPRECATED[upper]
            _warnings.warn(
                f"runlang='{self.runlang}' is deprecated. "
                f"Use '{new_equiv}' instead.",
                DeprecationWarning,
                stacklevel=3,
            )
            self._voter_snapshot_mode = snapshot
            return ""
        # Bare "C" defaults to C0 (final output)
        if upper == "C":
            self._voter_snapshot_mode = False
            return ""
        raise ValueError(
            f"Unknown Voter runlang code: '{self.runlang}'. "
            f"Valid: C0 (final), C0S (snapshots)."
        )

    def _build_c_arglist(self) -> list[str]:
        """Build argument list for unified VoterSimulator.

        The syshape argument carries the per-run dirname (Phase-C naming;
        see ``CBackendMixin._c_syshape_arg``), so every file the binary
        sprintf-composes — the transient ``s_*`` input and the
        ``m_*``/``sout_*``/``cldist_*`` observable outputs — lands inside
        ``self._c_rundir``; the read-back paths below use the exact same
        directory and C-side file names.
        """
        datdir = self._get_datdir_arg()
        syshape = self._c_syshape_arg()
        outdir = self._c_rundir if self._c_rundir is not None else self.dynpath
        self.out_id = self.out_suffix
        self.magn_path = outdir / self.sg.get_p_fname("m", self.out_id)

        # Track snapshot output path for cleanup
        if self._voter_snapshot_mode:
            self.sout_path = outdir / self.sg.get_p_fname(
                "sout",
                self.out_id,
            )
        else:
            self.sout_path = None

        # Track the cldist output path (sparse-binary; written by the C binary
        # when track_clusters, read back via _voter_io.load_cldist_bin).
        if self.track_clusters:
            self.cldist_path = outdir / self.sg.get_p_fname(
                "cldist",
                self.out_id,
            )
        else:
            self.cldist_path = None

        # Base (7) + rule/sampler/absorbing (6) + seed/clusters (3) = 16.
        arglist = [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.steps}",
            str(datdir),
            syshape,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_id),
            f"{RULE_CODE[self.rule]}",
            f"{self.q}",
            f"{self.eps:.12g}",
            f"{self.alpha:.12g}",
            f"{UPD_MODE_CODE[self.upd_mode]}",
            f"{1 if self.absorbing_check else 0}",
            f"{int(self.seed)}",
            f"{1 if self.track_clusters else 0}",
            f"{CLUSTER_MODE_CODE[self.cluster_mode]}",
        ]

        # 17th arg triggers snapshot mode in the unified binary
        if self._voter_snapshot_mode:
            arglist.append(f"{self.nSampleLog}")

        return arglist

    def run_cprogram(self, verbose: bool = False) -> None:
        """Execute C backend and read magnetization output."""
        # Call parent implementation for subprocess execution
        super().run_cprogram(verbose)
        # The C binary writes a magnetization series of length = sweeps actually
        # run (shorter than ``steps`` if it stopped early at an absorbing state).
        self.absorbed_at = None
        if self.magn_path and self.magn_path.exists():
            magn = np.fromfile(self.magn_path, dtype=np.float64)
            if self.absorbing_check and magn.size < self.steps:
                self.absorbed_at = int(magn.size) - 1
            # savemagn gates exposure (parity with py/pybind).
            if self.savemagn:
                self.observables[VOTER_OBS_MAGN].set_values(magn.tolist())
        # Cluster-size distribution (sparse-binary file from the C binary) ->
        # the same list[{size: count}] the py/pybind backends produce.
        if self.track_clusters:
            cldist_path = getattr(self, "cldist_path", None)
            if cldist_path is not None and Path(cldist_path).exists():
                self.observables[VOTER_OBS_CLDIST].set_values(
                    _voter_io.load_cldist_bin(cldist_path)
                )

    def _get_cleanup_paths(self) -> list[Path | None]:
        """Return paths to clean up after a C run.

        Under ``savedisk`` the C-written ``m_*``/``sout_*`` files are the
        persisted observables (stable ``out_suffix`` names), so they are kept;
        only the transient state/edgelist exports are removed.
        """
        paths: list[Path | None] = [getattr(self, "sfout", None)]
        if not self.savedisk:
            paths += [
                self.magn_path,
                getattr(self, "sout_path", None),
                getattr(self, "cldist_path", None),
            ]
        return paths

    def remove_run_c_files(self, remove_stderr: bool = True) -> None:
        """Remove the transient C exports; drop an empty non-savedisk rundir."""
        super().remove_run_c_files(remove_stderr)
        self._c_prune_rundir()

    # ------------------------------------------------------------------
    # Pybind11 backend (in-process, GT-compatible, seed-reproducible)
    # ------------------------------------------------------------------
    @staticmethod
    def _load_native_module():
        """Import the compiled ``_voter_native`` pybind11 module."""
        try:
            from .ccore import _voter_native  # type: ignore[import-untyped]

            return _voter_native
        except ImportError as exc:
            raise RuntimeError(
                "Pybind11 voter backend not available. Build the C extensions "
                "with `pip install -e .` or `make all` first."
            ) from exc

    def _run_pybind(self) -> None:
        """Run voter dynamics via the in-process pybind11 kernel.

        Reuses the same ``voter_apply_rule`` / ``voter_*_step`` C kernels as the
        subprocess backend (identical update logic across the rule family and
        sampler axis) but passes the graph as numpy CSR arrays, so there is no
        file I/O and GT graphs are supported.
        """
        mod = self._load_native_module()
        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        s_out, magn, absorbed_at, cluster_dist, s_t = mod.voter_sampling(
            self.s.astype(np.int8),
            ni,
            nw,
            nptr,
            int(self.steps),
            int(self.seed),
            bool(self.savemagn),
            int(RULE_CODE[self.rule]),
            int(self.q),
            float(self.eps),
            float(self.alpha),
            int(UPD_MODE_CODE[self.upd_mode]),
            bool(self.absorbing_check),
            bool(self.track_clusters),
            int(CLUSTER_MODE_CODE[self.cluster_mode]),
            bool(self.savedyn),
            (
                self._sout_sample_sweeps
                if self._sout_sample_sweeps is not None
                else np.empty(0, dtype=np.int64)
            ),
        )
        self.s = np.asarray(s_out, dtype=np.int8)
        self.observables[VOTER_OBS_MAGN].set_values(magn.tolist())
        self.absorbed_at = None if int(absorbed_at) < 0 else int(absorbed_at)
        if self.track_clusters:
            self.observables[VOTER_OBS_CLDIST].set_values(list(cluster_dist))
        if self.savedyn:
            # (n_recorded, N) int8 -> list of per-sweep 1-D arrays (matches py).
            # _persist_observables dumps this to sout_*.bin and drops the RAM
            # copy when savedisk=True (s_t then memmaps the file).
            self.observables[VOTER_OBS_SOUT].set_rows(
                [row.copy() for row in np.asarray(s_t, dtype=np.int8)]
            )

    # ------------------------------------------------------------------
    # Vectorized backends (NumPy `np` / CuPy `cu`)
    # ------------------------------------------------------------------
    def _assert_vectorized_supports_config(self, backend: str) -> None:
        """The vectorized backends implement the **synchronous** voter for the
        full rule family (linear / majority / qvoter / nonlinear) via per-sweep
        CSR gathers + segment reductions. They do not consult ``upd_mode``
        (always synchronous), and reject savedyn / the intrinsically-sequential
        schedules / cluster recording rather than silently mislabel a run.
        """
        if self.savedyn:
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not record the "
                f"per-sweep savedyn trajectory; use runlang='py'."
            )
        if self.upd_mode in ("link", "gillespie"):
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) is the synchronous "
                f"vectorized voter; upd_mode='{self.upd_mode}' is a different "
                f"schedule -- use runlang='py' or 'pb'."
            )
        if self.track_clusters:
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not emit the "
                f"cluster-size distribution (the fused vectorized sweep has no "
                f"per-sweep record hook) -- use runlang='py' or 'pb'."
            )

    def _run_vectorized(self, gpu: bool) -> None:
        """Run the vectorized synchronous voter (NumPy or CuPy).

        Every node updates from a frozen snapshot via per-sweep CSR gathers +
        segment reductions under the configured ``rule`` (linear / majority /
        qvoter / nonlinear); the schedule is synchronous and ``upd_mode`` is not
        consulted. Graph is passed as CSR arrays, so there is no file I/O and GT
        graphs are supported.
        """
        from ._vectorized_voter import run_vectorized_sync, CUPY_AVAILABLE

        if gpu:
            if not CUPY_AVAILABLE:
                raise RuntimeError(
                    f"CuPy backend requested (runlang='{self.runlang}') but cupy "
                    f"is not available; install cupy or use runlang='np'."
                )
            import cupy as xp
        else:
            xp = np
        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        s_out, magn, absorbed_at = run_vectorized_sync(
            xp,
            self.s.astype(np.int8),
            ni,
            nw,
            nptr,
            int(self.steps),
            int(self.seed),
            bool(self.savemagn),
            bool(self.absorbing_check),
            int(self.absorbing_every),
            rule=self.rule,
            q=int(self.q),
            eps=float(self.eps),
            alpha=float(self.alpha),
        )
        self.s = np.asarray(s_out, dtype=np.int8)
        self.observables[VOTER_OBS_MAGN].set_values(list(magn))
        self.absorbed_at = absorbed_at

    # ------------------------------------------------------------------
    # On-disk persistence of observables (savedisk; lazy disk-backed attrs)
    # ------------------------------------------------------------------
    def _is_py_runlang(self) -> bool:
        """True for the pure-Python loop (the only backend that streams sout)."""
        return not self.runlang.lower().startswith(("c", "pb", "np", "cu"))

    def _logspaced_sweeps(self, nlog: int) -> np.ndarray:
        """``nlog`` log-spaced sweep indices in ``[0, steps)`` (incl. 0, steps-1)."""
        steps = int(self.steps)
        if nlog >= steps:
            return np.arange(steps, dtype=np.int64)
        idx = np.unique(
            np.logspace(0.0, np.log10(steps - 1), int(nlog)).astype(np.int64)
        )
        idx = np.unique(np.concatenate(([0], idx, [steps - 1]))).astype(
            np.int64
        )
        return idx[(idx >= 0) & (idx < steps)]

    def _resolve_sout_sampling(self) -> "np.ndarray | None":
        """Choose the savedyn sampling: ``None`` = full (uniform ``sout_every``),
        otherwise a sorted array of log-spaced sweep indices.

        Falls back to ``VOTER_SOUT_NLOG_DEFAULT`` log-spaced snapshots (instead of
        raising) when the full trajectory would exceed ``VOTER_SOUT_MAX_BYTES``,
        unless ``sout_force_full`` overrides the cap. An explicit ``sout_nlog``
        always wins.
        """
        if not self.savedyn:
            return None
        steps = int(self.steps)
        n_full = -(-steps // self.sout_every)  # ceil(steps / sout_every)
        full_bytes = n_full * self.N  # int8 => 1 byte/elem
        if self.sout_nlog is not None and self.sout_nlog < steps:
            return self._logspaced_sweeps(self.sout_nlog)
        if full_bytes > VOTER_SOUT_MAX_BYTES:
            if self.sout_force_full:
                _warnings.warn(
                    f"sout_force_full=True: recording the full "
                    f"~{full_bytes / 2**20:.0f} MiB savedyn trajectory in memory "
                    f"(cap {VOTER_SOUT_MAX_BYTES / 2**20:.0f} MiB).",
                    stacklevel=3,
                )
                return None
            nlog = VOTER_SOUT_NLOG_DEFAULT
            _warnings.warn(
                f"savedyn full trajectory ~{full_bytes / 2**20:.0f} MiB exceeds "
                f"VOTER_SOUT_MAX_BYTES={VOTER_SOUT_MAX_BYTES / 2**20:.0f} MiB; "
                f"recording {nlog} log-spaced snapshots instead. Pass sout_nlog=K "
                f"to choose the count or sout_force_full=True to force the full "
                f"trajectory.",
                stacklevel=3,
            )
            return self._logspaced_sweeps(nlog)
        return None

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C: one directory per run)
    # ------------------------------------------------------------------
    def _lang_token(self) -> str:
        """Backend token for the run dirname.

        Resolved from ``runlang`` when ``run()`` has not stashed
        ``_active_backend`` yet: the manual ``init_voter_dynamics()`` flow
        exports the C inputs into the rundir BEFORE ``run()``, so the name
        must not shift between the export and the run.
        """
        backend = getattr(self, "_active_backend", None)
        if backend is None:
            backend = self._resolve_backend()
        return backend.value

    def _name_tokens(self) -> list:
        rule = self.rule
        return [
            ("p", float(self.sg.pflip)),
            ("q", self.q if rule == "qvoter" else None),
            ("eps", self.eps if rule in ("qvoter", "nonlinear") else None),
            ("alpha", self.alpha if rule == "nonlinear" else None),
            ("ns", int(self.steps)),
            ("rule", self.rule, DEFAULT_RULE),
            ("upd", self.upd_mode, DEFAULT_UPD_MODE),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _cfg_model_block(self) -> dict:
        """Every constructor argument with its resolved value."""
        return {
            "class": type(self).__name__,
            "rule": self.rule,
            "q": self.q,
            "eps": self.eps,
            "alpha": self.alpha,
            "upd_mode": self.upd_mode,
            "absorbing_check": self.absorbing_check,
            "absorbing_every": self.absorbing_every,
            "track_clusters": self.track_clusters,
            "cluster_mode": self.cluster_mode,
            "freq": self.freq,
            "nSampleLog": self.nSampleLog,
            "sout_every": self.sout_every,
            "sout_nlog": self.sout_nlog,
            "sout_force_full": self.sout_force_full,
            "savemagn": self.savemagn,
            "savedyn": self.savedyn,
            "savedisk": self.savedisk,
            "ic": self.ic,
            "state_type": self.state_type,
        }

    def _begin_outputs(self) -> None:
        """Resolve the per-run output directory and open the sout stream.

        Layout (Phase C): ONE directory per run under ``dynpath``, named by
        the canonical parameter tokens (``_name_tokens``); the in-process
        backends write the SHORT canonical names (``m.bin`` / ``sout.bin``
        / ``cldist.npz``) plus the ``cfg.json`` sidecar inside it, while the
        C subprocess writes its own sprintf-composed names into the SAME
        rundir (redirected via the composed syshape argument). No-op when
        ``savedisk=False`` or when nothing is persisted (lazy filesystem).
        """
        sout = self.observables[VOTER_OBS_SOUT]
        # Always start from a clean streaming state.
        sout.reset_stream()
        self._sout_seen = 0
        # Resolve savedyn sampling (full | log-spaced fallback | forced-full) for
        # the in-process backends, regardless of savedisk; the C subprocess samples
        # its own file snapshots via nSampleLog.
        self._sout_sample_sweeps = (
            self._resolve_sout_sampling()
            if self.savedyn and not self.runlang.startswith("C")
            else None
        )
        self._sout_sample_set = (
            set(self._sout_sample_sweeps.tolist())
            if self._sout_sample_sweeps is not None
            else None
        )
        if not self.savedisk:
            return
        if self.runlang.startswith("C"):
            # The C exchange targets the rundir stashed when the arglist was
            # built (identical tokens); read-back paths were set there too.
            # Only the reproducibility sidecar is written from Python.
            rundir = (
                self._c_rundir
                if self._c_rundir is not None
                else self._run_output_dir()
            )
            rundir.mkdir(parents=True, exist_ok=True)
            self._write_cfg_sidecar(rundir)
            return
        if not (self.savemagn or self.track_clusters or self.savedyn):
            return  # nothing persisted -> no directory (lazy filesystem)
        rundir = self._run_output_dir()
        rundir.mkdir(parents=True, exist_ok=True)
        self._write_cfg_sidecar(rundir)
        if self.savemagn:
            self.observables[VOTER_OBS_MAGN].set_path(
                rundir / f"{VOTER_MAGN_FBASE}{BIN}"
            )
        if self.track_clusters:
            self.observables[VOTER_OBS_CLDIST].set_path(
                rundir / f"{VOTER_CLDIST_FBASE}{NPZ}"
            )
        if self.savedyn:
            sout.set_path(rundir / f"{VOTER_SOUT_FBASE}{BIN}")
            if self._is_py_runlang():
                sout.open_stream()

    def _persist_observables(self) -> None:
        """Finalize disk artifacts and switch attributes to disk-backed views."""
        if not self.savedisk:
            return
        magn = self.observables[VOTER_OBS_MAGN]
        sout = self.observables[VOTER_OBS_SOUT]
        cldist = self.observables[VOTER_OBS_CLDIST]
        # Close a streamed sout (Python loop) -> finalize sout_*.bin (no-op if the
        # stream was never opened, e.g. the pybind / C backends).
        sout.close_stream()
        if self.runlang.startswith("C"):
            # The C binary already wrote m_*/sout_* under the stable names; just
            # point the lazy views at them (cleanup is suppressed via
            # _get_cleanup_paths). C does not emit cldist.
            if self.savemagn:
                magn.set_path(self.magn_path)
            sout_p = getattr(
                self, "sout_path", None
            )  # set in C0S snapshot mode
            if sout_p is not None and Path(sout_p).exists():
                sout.set_path(sout_p)
                sout.mark_on_disk()
            return
        # In-process backends: dump the RAM buffers. For pybind the whole
        # (n_rec, N) array is cached (dumped, subsampled, then the RAM copy
        # dropped); the Python loop already streamed its rows (cache empty ->
        # finalize_batch only drops it once the file exists).
        if self.savedyn and sout.path is not None:
            sout.finalize_batch(subsample=self.sout_every)
        if self.savemagn:
            magn.persist()
        if self.track_clusters:
            cldist.persist()

    # -- public loaders / introspection (read big output only when needed) --
    def _require_output(self, name: str, what: str) -> Path:
        path = self.observables[name].path
        if path is None or not Path(path).exists():
            raise FileNotFoundError(
                f"No {what} file on disk. Run with savedisk=True and the matching "
                f"save flag set first (and check out_suffix)."
            )
        return Path(path)

    def load_magn(self) -> np.ndarray:
        """Load the persisted magnetization series (``float64`` array)."""
        return _voter_io.load_magn(
            self._require_output(VOTER_OBS_MAGN, "magnetization")
        )

    def load_sout(self, mmap: bool = False) -> np.ndarray:
        """Load the persisted spin trajectory as ``(n_rec, N)`` int8.

        ``mmap=True`` returns a read-only memmap (no full-RAM copy).
        """
        return _voter_io.load_sout(
            self._require_output(VOTER_OBS_SOUT, "sout trajectory"),
            self.N,
            mmap=mmap,
        )

    def load_clusters(self) -> list[dict[int, int]]:
        """Load the persisted cluster-size-distribution time series."""
        return _voter_io.load_cldist(
            self._require_output(VOTER_OBS_CLDIST, "cluster-size distribution")
        )

    def output_sizes(self) -> dict[str, int]:
        """Bytes on disk for each persisted observable (0 if not written)."""
        sizes: dict[str, int] = {}
        for name in (VOTER_OBS_MAGN, VOTER_OBS_SOUT, VOTER_OBS_CLDIST):
            path = self.observables[name].path
            sizes[name] = (
                Path(path).stat().st_size
                if path is not None and Path(path).exists()
                else 0
            )
        return sizes

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def run(
        self,
        tqdm_on: bool = False,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        """Run voter model dynamics.

        Thin wrapper over the shared ``RunHostMixin.run`` skeleton (this
        class's hand-copied version was hoisted into the mixin); it only
        resolves the legacy ``eqSTEP`` alias (``steps`` wins when both are
        given, matching ``initialize_run_parameters``).

        Parameters
        ----------
        tqdm_on : bool
            Show progress bar.
        steps : int, optional
            Override number of steps.
        simref : float, optional
            Size-normalized time.
        eqSTEP : int, optional
            Legacy alias for steps.
        verbose : bool
            Verbose output.
        clean_export : bool
            Remove the C backend's transient export files (state/edgelist) after
            a C run. Persisted observable files (``m_*``/``sout_*``) are kept when
            ``savedisk=True``; only suppressed in-RAM mode removes them.

        Notes
        -----
        With ``savedisk=True`` (default), enabled observables are written under
        ``dynpath`` and ``magn`` / ``s_t`` / ``cluster_dist`` become lazy,
        disk-backed views; see the class docstring and ``load_sout`` /
        ``output_sizes``.
        """
        super().run(
            tqdm_on=tqdm_on,
            steps=steps if steps is not None else eqSTEP,
            simref=simref,
            verbose=verbose,
            clean_export=clean_export,
        )


# Register VoterModel's solver backends (py/pb/np/cu/C) in the shared solver
# registry. Imported after the class so that importing VoterModel wires up
# dispatch; no circular import (``_solvers`` does not import this class).
from . import _solvers  # noqa: E402,F401
