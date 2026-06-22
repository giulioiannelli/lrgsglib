"""Voter-model dynamics on signed graphs.

Runlang convention: ``C<digit><letters>``
  - digit 0 = voter dynamics (single algorithm)
  - letters: S = snapshots (bare digit = final state only)
"""

from __future__ import annotations

import warnings as _warnings
from pathlib import Path
from typing import Any

import numpy as np
import tqdm

from .._c_backend import CBackendMixin
from .._csr import build_graph_csr
from ..BinDynSys import BinDynSys
from ...utils.tools.chronometer import time_function_accumulate
from .defaults import (
    CLUSTER_MODE_CODE,
    CLUSTER_MODES,
    CLUSTER_SPLIT_SCAN_CAP,
    DEFAULT_ABSORBING_EVERY,
    DEFAULT_CLUSTER_MODE,
    DEFAULT_NOISE_EPS,
    DEFAULT_NONLIN_ALPHA,
    DEFAULT_QVOTER_Q,
    DEFAULT_RULE,
    DEFAULT_UPD_MODE,
    GILLESPIE_RULES,
    LINK_ONLY_RULES,
    PY_ONLY_UPD_MODES,
    RULE_CODE,
    UPD_MODE_CODE,
    VOTER_RULES,
    VOTER_UPD_MODES,
    VOTER_UPD_MODES_PLANNED,
)

# Update schedules for which the per-flip cluster-size-distribution tracker is
# wired today (single-spin engines). synchronous / link / vectorized are LATER.
_CLUSTER_UPD_MODES: frozenset[str] = frozenset({"asynchronous", "gillespie"})

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


class VoterModel(CBackendMixin, BinDynSys):
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
    save_magnetization : bool, optional
        If True, record the magnetization series ``magn`` at each sweep.
        Honoured identically by the Python and C backends.
    upd_mode : str, optional
        Update schedule (Axis B): ``'asynchronous'`` (random sequential),
        ``'synchronous'`` (double-buffer), ``'link'`` (edge-update, linear
        only), or ``'gillespie'`` (rejection-free CTMC, linear only). All
        schedules are available on every backend (Python, C subprocess, pybind).
    freq : int, optional
        Recording frequency.
    nSampleLog : int, optional
        Number of log-spaced samples (for C backend).
    **kwargs
        Additional arguments passed to BinDynSys.

    Examples
    --------
    >>> voter = VoterModel(lattice, steps=100, runlang='py')
    >>> voter.init_voter_dynamics()
    >>> voter.run(tqdm_on=False)
    """

    dyn_UVclass = "voter_model"

    # CBackendMixin configuration
    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = "VoterSimulator{}"
    # Validation handled by suffix method
    _allowed_c_keys: tuple[str, ...] = ()

    # Class-level observable defaults (will be overwritten per instance)
    magn: list[float] = []
    s_t: list[np.ndarray] = []
    _voter_snapshot_mode: bool = False

    def __init__(
        self,
        sg,
        *,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
        save_magnetization: bool = False,
        rule: str = DEFAULT_RULE,
        q: int = DEFAULT_QVOTER_Q,
        eps: float = DEFAULT_NOISE_EPS,
        alpha: float = DEFAULT_NONLIN_ALPHA,
        upd_mode: str = DEFAULT_UPD_MODE,
        absorbing_check: bool = False,
        absorbing_every: int = DEFAULT_ABSORBING_EVERY,
        track_clusters: bool = False,
        cluster_mode: str = DEFAULT_CLUSTER_MODE,
        freq: int = 10,
        nSampleLog: int = 100,
        **kwargs: Any,
    ) -> None:
        dynpath = getattr(sg, 'path_voter', None)
        resolved_steps = steps if steps is not None else eqSTEP
        super().__init__(sg, dynpath=dynpath, steps=resolved_steps, simref=simref, **kwargs)
        self.save_magnetization = save_magnetization
        self.rule = self._validate_rule(rule)
        self.q, self.eps, self.alpha = self._validate_rule_params(q, eps, alpha)
        self.upd_mode = self._validate_upd_mode(upd_mode)
        self._validate_rule_mode(self.rule, self.upd_mode)
        self.absorbing_check = bool(absorbing_check)
        self.absorbing_every = max(1, int(absorbing_every))
        self.absorbed_at: int | None = None
        self.track_clusters = bool(track_clusters)
        self.cluster_mode = self._validate_cluster_mode(cluster_mode)
        self._tracker = None
        self.freq = freq
        self.nSampleLog = nSampleLog
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
        raise ValueError(
            f"Unknown rule='{rule}'. Valid: {list(VOTER_RULES)}."
        )

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
        q: int, eps: float, alpha: float,
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
        """Reset cached observables collected during a run."""
        self.magn = []
        self.s_t = []
        self.cluster_dist: list[dict[int, int]] = []
        self.absorbed_at = None

    def init_voter_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the spin configuration and export data if required."""
        self._check_c_backend_or_fallback()
        self.reset_observables()
        self.init_s(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            # Export uses build_p_fname which handles underscore via join_non_empty
            self.sg._export_edgel_bin(exName=self.run_id)
            self.sg.export_adj_bin(exName=self.run_id)
        self.sini = self.s.copy()

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
            return np.int8(1 if np.random.randint(2) else -1)  # random tie-break
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
            # majority. Project-chosen form (citation pending; see defaults.py).
            deg = len(nbrs)
            nplus = 0
            for j, w in nbrs:
                if ((-1 if w < 0 else 1) * int(s_read[j])) == 1:
                    nplus += 1
            fp = nplus / deg
            fm = 1.0 - fp
            num = fp ** self.alpha
            den = num + fm ** self.alpha
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
        return sum(1 for u, v, w in self._edges() if w * int(s[u]) * int(s[v]) < 0)

    def is_absorbing(self, s: np.ndarray | None = None) -> bool:
        """True iff ``s`` has zero frustrated edges (the dynamics is frozen).

        A frozen voter configuration exists iff the signed graph is balanced
        (``lambda_min(signed Laplacian) ~ 0``); on a frustrated graph it never
        occurs and the dynamics never freezes (ref [3], Properties (ii)-(iii),
        p. 9). On an unsigned graph this reduces to consensus ``|M| = 1``.
        """
        return self.count_frustrated_edges(s) == 0

    # -- per-sweep samplers (Axis B) ----------------------------------
    def _record(self) -> None:
        if self.save_magnetization:
            self.magn.append(float(np.sum(self.s)) / float(self.N))
        if self.savedyn:
            self.s_t.append(self.s.copy())
        if self.track_clusters and self._tracker is not None:
            self.cluster_dist.append(self._tracker.size_distribution())

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
        tr = self._tracker
        if tr is None:
            for _ in range(N):
                nd = int(np.random.randint(N))
                s[nd] = self._apply_rule(nd, s)
        else:
            # Cluster tracking: feed every genuine sign flip to the tracker
            # (a non-flipping update leaves the partition untouched).
            for _ in range(N):
                nd = int(np.random.randint(N))
                new = self._apply_rule(nd, s)
                if new != s[nd]:
                    tr.flip(nd)        # pre-flip bookkeeping
                    s[nd] = new

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

    # -- cluster-size-distribution tracker (single-spin engines) ------
    def _setup_cluster_tracker(self) -> None:
        """Build the incremental cluster tracker for the Python single-spin
        schedules (asynchronous / gillespie), or clear it when not tracking.

        synchronous / link (and the vectorized backends) are not wired yet and
        raise rather than silently dropping the observable.
        """
        if not self.track_clusters:
            self._tracker = None
            return
        if self.upd_mode not in _CLUSTER_UPD_MODES:
            raise NotImplementedError(
                f"track_clusters is wired only for upd_mode in "
                f"{sorted(_CLUSTER_UPD_MODES)} (single-spin engines); "
                f"got upd_mode='{self.upd_mode}'."
            )
        from ._cluster_tracker import ClusterTracker, edge_sign_arrays
        idx, signs = self._gillespie_neighbors()
        b = edge_sign_arrays(signs, self.cluster_mode)
        self._tracker = ClusterTracker(self.s, idx, b)

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
        t = 0.0          # continuous time, in sweep units
        recorded = 0     # next integer sweep-time still to record
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
                    while recorded < steps:   # pad the frozen tail
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
            i = int(np.searchsorted(np.cumsum(rate), np.random.random() * R,
                                    side="right"))
            if i >= N:
                i = N - 1
            if self._tracker is not None:
                self._tracker.flip(i)        # pre-flip bookkeeping
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
        self._setup_cluster_tracker()
        if self.upd_mode == "gillespie":
            self._run_gillespie(tqdm_on)
            return
        sweep_fn = {
            "asynchronous": self._sweep_async,
            "synchronous": self._sweep_sync,
            "link": self._sweep_link,
        }[self.upd_mode]
        iterator = tqdm.tqdm(range(self.steps)) if tqdm_on else range(self.steps)
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
        absorbing early-stop, but capture only the final state and the
        magnetization series -- not the per-sweep ``savedyn`` trajectory.
        Refuse rather than silently return an empty ``s_t``. ``PY_ONLY_UPD_MODES``
        (currently empty) lists any schedule that is still Python-only.
        """
        if self.savedyn:
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not record the "
                f"per-sweep savedyn trajectory; use runlang='py' (full per-sweep "
                f"s_t) or 'C0S' (nSampleLog file snapshots)."
            )
        if self.upd_mode in PY_ONLY_UPD_MODES:
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not implement "
                f"upd_mode='{self.upd_mode}' yet (native CTMC kernel pending); "
                f"use runlang='py'."
            )
        if self.track_clusters and not (
            backend == "pybind" and self.upd_mode == "gillespie"
        ):
            # Native cluster tracking is wired for pybind + gillespie (the shared
            # _ccore CTMC kernel); the C subprocess does not emit it yet.
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) emits the cluster-size "
                f"distribution only via runlang='pb_voter' with "
                f"upd_mode='gillespie'; use runlang='py' otherwise."
            )

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
                DeprecationWarning, stacklevel=3,
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
        """Build argument list for unified VoterSimulator."""
        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except ValueError:
            datdir = self.sg.path_sgdata
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        self.out_id = self.out_suffix
        self.magn_path = self.dynpath / self.sg.get_p_fname('m', self.out_id)

        # Track snapshot output path for cleanup
        if self._voter_snapshot_mode:
            self.sout_path = self.dynpath / self.sg.get_p_fname(
                'sout', self.out_id,
            )
        else:
            self.sout_path = None

        # Base arguments (7) + Phase-3 rule/sampler/absorbing (6) = 13.
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
        ]

        # 14th arg triggers snapshot mode in the unified binary
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
            # save_magnetization gates exposure (parity with py/pybind).
            if self.save_magnetization:
                self.magn = magn.tolist()

    def _get_cleanup_paths(self) -> list[Path | None]:
        """Return paths to clean up after C run."""
        return [
            getattr(self, 'sfout', None),
            self.magn_path,
            getattr(self, 'sout_path', None),
        ]

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
        s_out, magn, absorbed_at, cluster_dist = mod.voter_sampling(
            self.s.astype(np.int8),
            ni, nw, nptr,
            int(self.steps),
            int(self.seed),
            bool(self.save_magnetization),
            int(RULE_CODE[self.rule]),
            int(self.q),
            float(self.eps),
            float(self.alpha),
            int(UPD_MODE_CODE[self.upd_mode]),
            bool(self.absorbing_check),
            bool(self.track_clusters),
            int(CLUSTER_MODE_CODE[self.cluster_mode]),
            int(CLUSTER_SPLIT_SCAN_CAP),
        )
        self.s = np.asarray(s_out, dtype=np.int8)
        self.magn = magn.tolist()
        self.absorbed_at = None if int(absorbed_at) < 0 else int(absorbed_at)
        if self.track_clusters:
            self.cluster_dist = list(cluster_dist)

    # ------------------------------------------------------------------
    # Vectorized backends (NumPy `np_voter` / CuPy `cu_voter`)
    # ------------------------------------------------------------------
    def _assert_vectorized_supports_config(self, backend: str) -> None:
        """The vectorized backends implement the **synchronous linear** voter
        only -- a single CSR gather per sweep. They do not consult ``upd_mode``
        (always synchronous), and reject the rule family / savedyn / the
        intrinsically-sequential schedules rather than silently mislabel a run.
        """
        if self.rule != "linear":
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) implements the vectorized "
                f"synchronous LINEAR voter only; rule='{self.rule}' is "
                f"unsupported -- use runlang in (py, C0*, pb_voter)."
            )
        if self.savedyn:
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not record the "
                f"per-sweep savedyn trajectory; use runlang='py'."
            )
        if self.upd_mode in ("link", "gillespie"):
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) is the synchronous "
                f"vectorized voter; upd_mode='{self.upd_mode}' is a different "
                f"schedule -- use runlang='py' or 'pb_voter'."
            )
        if self.track_clusters:
            raise NotImplementedError(
                f"runlang='{self.runlang}' ({backend}) does not emit the "
                f"cluster-size distribution (synchronous vectorized; tracker is "
                f"single-spin) -- use runlang='py'."
            )

    def _run_vectorized(self, gpu: bool) -> None:
        """Run the vectorized synchronous linear voter (NumPy or CuPy).

        Every node copies a uniformly chosen signed neighbour from a frozen
        snapshot via one CSR gather per sweep (the schedule is synchronous;
        ``upd_mode`` is not consulted). Graph is passed as CSR arrays, so there
        is no file I/O and GT graphs are supported.
        """
        from ._vectorized_voter import run_vectorized_sync, CUPY_AVAILABLE
        if gpu:
            if not CUPY_AVAILABLE:
                raise RuntimeError(
                    f"CuPy backend requested (runlang='{self.runlang}') but cupy "
                    f"is not available; install cupy or use runlang='np_voter'."
                )
            import cupy as xp
        else:
            xp = np
        ni, nw, nptr = build_graph_csr(self.sg, self.N)
        s_out, magn, absorbed_at = run_vectorized_sync(
            xp, self.s.astype(np.int8), ni, nw, nptr,
            int(self.steps), int(self.seed),
            bool(self.save_magnetization),
            bool(self.absorbing_check), int(self.absorbing_every),
        )
        self.s = np.asarray(s_out, dtype=np.int8)
        self.magn = magn
        self.absorbed_at = absorbed_at

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        """Run voter model dynamics.

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
            Remove exported files after run.
        """
        is_c = self.runlang.startswith("C")
        is_pb = self.runlang.lower().startswith("pb")
        is_np = self.runlang.lower().startswith("np")
        is_cu = self.runlang.lower().startswith("cu")
        if is_c or is_pb:
            # Guard before any file export / kernel call so the native backends
            # never silently ignore rule / upd_mode / absorbing_check.
            self._assert_native_supports_config("C subprocess" if is_c else "pybind")
        if is_np or is_cu:
            self._assert_vectorized_supports_config("CuPy" if is_cu else "NumPy vectorized")
        self.check_attribute()
        self.initialize_run_parameters(steps=steps, simref=simref, eqSTEP=eqSTEP)
        if is_c:
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        elif is_pb:
            self._run_pybind()
        elif is_np:
            self._run_vectorized(gpu=False)
        elif is_cu:
            self._run_vectorized(gpu=True)
        else:
            self.voter_sampling(tqdm_on)
