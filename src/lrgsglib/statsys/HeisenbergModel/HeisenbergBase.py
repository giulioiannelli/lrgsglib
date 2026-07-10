"""HeisenbergBase — the classical O(3) substrate on the shared statsys spine.

The Layer-0 model class of the Heisenberg modernization (statsys unification
D-B8/D-B9; the second continuous-state consumer after XY, now with a VECTOR
state — one unit 3-vector per site): it owns the Hamiltonian (ONE symmetric
coupling matrix ``J`` built once by the shared
``statsys._csr.build_coupling_csr``, from which BOTH ``delta_E`` and the
energy observable derive — ΔE↔E consistency by construction, plan §3b), the
substrate contracts consumed by the shared ``statsys._thermal`` engines
(single-site + cluster + exchange), and the ObservableSet lifecycle.

It owns NOTHING scheme-specific — temperature policy, acceptance rule and
update axes live on the thin scheme leaves (``HeisenbergMetropolis``, ...).

This class is NEW code beside the legacy ``HeisenbergModel`` (strangler
pattern): the legacy class keeps working untouched; both dispatch through the
same ``"heisenberg"`` solver-registry key, and the solver polymorphs on the
instance.

Hamiltonian (classical O(3); ``k_B = 1``)::

    H = − Σ_{(i,j) ∈ E} J_ij (n_i · n_j),      |n_i| = 1
    ΔE(i: n → n') = (n − n') · Σ_j J_ij n_j

with ``J`` from the signed edge weights via ``coupling_norm``. The external-
field term ``−h Σ (n_i · e_h)`` is NOT wired: the field must be zero (hard
capability error, invariant #3).

Cluster moves use the embedded-Ising construction of ref M6 (Wolff 1989, the
O(n) case): draw a random unit vector ``r``, grow FK clusters on the bond
satisfaction ``sat_ij = J_ij (n_i · r)(n_j · r)`` and flip a cluster by
reflecting its spins about the hyperplane perpendicular to ``r``,
``n → n − 2(n·r) r``. A reflection is an involution and flips the sign of
``n·r``, so the satisfaction is flip-covariant exactly as in the Ising case —
cluster moves therefore stay valid on SIGNED couplings.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np
import tqdm

from ...config.const import BIN
from .._c_backend import CBackendMixin
from .._csr import build_coupling_csr
from .._observables import ObservableSet, RowTrajectory, ScalarSeries
from .._run_host import RunHostMixin
from .._vectorized import VectorSyncEngine  # noqa: F401  (leaves build it)
from ..VecDynSys import VecDynSys
from .defaults import (
    HEISENBERG_COUPLING_NORM_DEFAULT,
    HEISENBERG_COUPLING_NORMS,
    HEISENBERG_DELTA_DEFAULT,
    HEISENBERG_DYN_SUBDIR,
    HEISENBERG_ENERGY_FBASE,
    HEISENBERG_MAGN_FBASE,
    HEISENBERG_OBS_DEFAULT,
    HEISENBERG_OBS_ENERGY,
    HEISENBERG_OBS_MAGN,
    HEISENBERG_OBS_SNAPSHOTS,
    HEISENBERG_SNAPSHOT_EVERY_DEFAULT,
    HEISENBERG_SNAPSHOT_MAX_BYTES,
    HEISENBERG_SNAPSHOTS_FBASE,
    HEISENBERG_SOLVER_NAME,
    HEISENBERG_T_DEFAULT,
)

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["HeisenbergBase"]

_KNOWN_OBSERVABLES = (
    HEISENBERG_OBS_ENERGY,
    HEISENBERG_OBS_MAGN,
    HEISENBERG_OBS_SNAPSHOTS,
)


class HeisenbergBase(RunHostMixin, CBackendMixin, VecDynSys):
    """O(3) substrate: symmetric-J Hamiltonian + observables + run host.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings ``w_ij``.
    T : float
        Temperature (``k_B = 1``; stored by ``VecDynSys``, the policy lives
        on the scheme leaves).
    delta : float
        Perturbation amplitude of the single-site proposal:
        ``normalize(n + delta·u)`` with ``u`` uniform in [−1, 1]³ — the
        historical kernel shared with the legacy and native backends.
    coupling_norm : {'raw', 'sym', 'avg'}
        How the symmetric ``J`` derives from the edge weights (shared
        vocabulary, ``statsys._csr``).
    observables : sequence of str
        Which series to record: any of ``'energy'`` (intensive per-site
        energy, tracked incrementally from the sweep ΔE), ``'magn'`` (the
        magnetisation ``|Σ n_i|/N``) and ``'snapshots'`` ((n_rec, 3N)
        float64 trajectory — each row is the flattened (N, 3) state,
        subsampled by *snapshot_every*).
    snapshot_every : int
        Keep every k-th sweep in the snapshot trajectory.
    **kwargs
        Forwarded to :class:`VecDynSys` (``ic``, ``runlang``,
        ``steps``/``simref``, ``seed``, ``savedisk``, ...).
    """

    solver_name = HEISENBERG_SOLVER_NAME

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        T: float = HEISENBERG_T_DEFAULT,
        delta: float = HEISENBERG_DELTA_DEFAULT,
        coupling_norm: str = HEISENBERG_COUPLING_NORM_DEFAULT,
        observables: Sequence[str] = HEISENBERG_OBS_DEFAULT,
        snapshot_every: int = HEISENBERG_SNAPSHOT_EVERY_DEFAULT,
        **kwargs: Any,
    ) -> None:
        if not (float(delta) > 0.0):
            raise ValueError(f"delta must be > 0; got {delta}.")
        if coupling_norm not in HEISENBERG_COUPLING_NORMS:
            raise ValueError(
                f"coupling_norm must be one of {HEISENBERG_COUPLING_NORMS}; "
                f"got {coupling_norm!r}."
            )
        selected = tuple(observables)
        unknown = [o for o in selected if o not in _KNOWN_OBSERVABLES]
        if unknown:
            raise ValueError(
                f"Unknown observable(s) {unknown}; known: {_KNOWN_OBSERVABLES}."
            )
        self.observables = ObservableSet(
            [
                ScalarSeries(HEISENBERG_OBS_ENERGY, dtype=np.float64),
                ScalarSeries(HEISENBERG_OBS_MAGN, dtype=np.float64),
                # One row per record = the flattened (N, 3) state.
                RowTrajectory(
                    HEISENBERG_OBS_SNAPSHOTS, ncols=3 * sg.N, dtype=np.float64
                ),
            ]
        )
        self._selected_obs = selected
        self.delta = float(delta)
        self.coupling_norm = coupling_norm
        self.snapshot_every = max(1, int(snapshot_every))

        # Graph-anchored dynamics tree (data/<graph>/<model>/N=...),
        # like path_ising; legacy classes keep the old flat location.
        dynpath = getattr(sg, "path_heisenberg", None)
        if dynpath is None:
            base = getattr(sg, "path_data", None)
            if base is not None:
                dynpath = Path(base) / HEISENBERG_DYN_SUBDIR
        super().__init__(
            sg,
            q=3,
            discrete=False,
            T=float(T),
            dynpath=dynpath,
            **kwargs,
        )
        # The Heisenberg field term −h Σ (n_i · e_h) needs a field DIRECTION
        # on top of DynSys's per-site amplitude vector; it is not wired.
        if np.any(np.asarray(self.field, dtype=np.float64) != 0.0):
            raise NotImplementedError(
                "The Heisenberg external-field term −h Σ (n_i · e_h) is not "
                "wired; the field must be zero."
            )

        # Configure-then-run: the recorder list is compiled once here; the
        # per-sweep _record only iterates it.
        recorders = []
        if HEISENBERG_OBS_ENERGY in selected:
            recorders.append(self._record_energy)
        if HEISENBERG_OBS_MAGN in selected:
            recorders.append(self._record_magn)
        if HEISENBERG_OBS_SNAPSHOTS in selected:
            recorders.append(self._record_snapshot)
        self._recorders = tuple(recorders)

        # Substrate arrays (built once in init_dynamics)
        self._nbr_idx: np.ndarray | None = None
        self._nbr_J: np.ndarray | None = None
        self._nbr_ptr: np.ndarray | None = None
        self._nbr_rows: np.ndarray | None = None
        self._edge_u: np.ndarray | None = None
        self._edge_v: np.ndarray | None = None
        self._edge_J: np.ndarray | None = None
        # Embedded-Ising reflection direction for the cluster moves (M6):
        # the leaf draws it per cluster (wolff) or per sweep (sw).
        self._refl_dir: np.ndarray = np.array([0.0, 0.0, 1.0])
        self._redraw_dir_on_flip: bool = False
        # Running extensive energy, advanced by the engine's per-sweep ΔE.
        self._e_running: float | None = None
        # The engine of the last run (per-move extras live on it).
        self._engine: Any = None
        self.sini: np.ndarray | None = None
        self.reset_observables()

    @property
    def _state_shape(self) -> tuple[int, ...]:
        return (self.N, 3)

    # ------------------------------------------------------------------
    # State initialisation (unit 3-vectors on S²)
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise unit vectors: 'uniform' draws isotropically on the
        sphere (normalized Gaussians); 'zero' aligns everything along +z;
        'custom' takes an (N, 3) array and normalizes its rows."""
        match self.ic:
            case "uniform" | "random" | "rand":
                vecs = np.random.randn(self.N, 3)
                norms = np.linalg.norm(vecs, axis=1, keepdims=True)
                norms = np.where(norms > 0, norms, 1.0)
                self.s = (vecs / norms).astype(np.float64)
            case "zero" | "zeros":
                self.s = np.zeros((self.N, 3), dtype=np.float64)
                self.s[:, 2] = 1.0
            case "custom":
                if custom is None:
                    raise ValueError("Provide a custom state array.")
                self.s = np.asarray(custom, dtype=np.float64).copy()
                if self.s.shape != (self.N, 3):
                    raise ValueError(
                        f"Shape mismatch: {self.s.shape} != ({self.N}, 3)"
                    )
                norms = np.linalg.norm(self.s, axis=1, keepdims=True)
                self.s /= np.where(norms > 0, norms, 1.0)
            case _:
                raise ValueError(f"Unknown IC: {self.ic!r}")

    # ------------------------------------------------------------------
    # Disk-backed observables (lazy views over the ObservableSet)
    # ------------------------------------------------------------------
    @property
    def ene(self) -> list:
        """Intensive per-site energy series e(t) = E(t)/N."""
        return self.observables[HEISENBERG_OBS_ENERGY].data

    @ene.setter
    def ene(self, value) -> None:
        self.observables[HEISENBERG_OBS_ENERGY].set_values(value)

    @property
    def magn(self) -> list:
        """Magnetisation series m(t) = |Σ n_i|/N."""
        return self.observables[HEISENBERG_OBS_MAGN].data

    @magn.setter
    def magn(self, value) -> None:
        self.observables[HEISENBERG_OBS_MAGN].set_values(value)

    @property
    def s_t(self):
        """Spin trajectory (rows = flattened (N, 3) states): RAM list
        pre-persist, read-only memmap after."""
        return self.observables[HEISENBERG_OBS_SNAPSHOTS].data

    @s_t.setter
    def s_t(self, value) -> None:
        self.observables[HEISENBERG_OBS_SNAPSHOTS].set_rows(value)

    # ------------------------------------------------------------------
    # Couplings (ONE symmetric J, built once — plan §3b, shared builder)
    # ------------------------------------------------------------------
    def _ensure_couplings(self) -> None:
        if self._nbr_J is not None:
            return
        c = build_coupling_csr(self.sg, self.N, self.coupling_norm)
        self._nbr_idx = c.idx
        self._nbr_J = c.J
        self._nbr_ptr = c.ptr
        self._nbr_rows = c.rows
        # Each undirected edge once (u < v), for the SW bond sweep + energy.
        self._edge_u = c.edge_u
        self._edge_v = c.edge_v
        self._edge_J = c.edge_J

    # ------------------------------------------------------------------
    # Substrate contract (consumed by statsys._thermal.ThermalEngine)
    # ------------------------------------------------------------------
    def propose_local(self, nd: int) -> np.ndarray:
        """``normalize(n + delta·u)`` with ``u`` uniform in [−1, 1]³ — the
        historical proposal shared with the legacy loop and the native
        kernel."""
        perturb = self.s[nd] + self.delta * np.random.uniform(-1.0, 1.0, 3)
        norm = float(np.linalg.norm(perturb))
        if norm > 0.0:
            return perturb / norm
        return self.s[nd].copy()

    def _local_field(self, nd: int) -> np.ndarray:
        """``Σ_j J_nd,j n_j`` — the molecular field at site *nd*."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        return self._nbr_J[lo:hi] @ self.s[self._nbr_idx[lo:hi]]

    def delta_E(self, nd: int, proposal: np.ndarray) -> float:
        """``(n − n') · Σ_j J_nd,j n_j``."""
        return float(np.dot(self.s[nd] - proposal, self._local_field(nd)))

    def commit(self, nd: int, proposal: np.ndarray) -> None:
        self.s[nd] = proposal

    def order_parameter(self, s: np.ndarray | None = None) -> float:
        """Magnetisation ``m = |Σ n_i|/N`` — 1 for aligned spins, of order
        ``N^{−1/2}`` in the disordered phase."""
        state = self.s if s is None else np.asarray(s)
        return float(np.linalg.norm(np.mean(state, axis=0)))

    # --- cluster contract (statsys._thermal.ClusterSubstrate) ---
    def _draw_reflection(self) -> None:
        """Draw a fresh embedded-Ising reflection direction r (ref M6):
        isotropic on the unit sphere."""
        while True:
            r = np.random.randn(3)
            norm = float(np.linalg.norm(r))
            if norm > 1e-12:
                self._refl_dir = r / norm
                return

    def neighbor_indices(self, nd: int) -> np.ndarray:
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        return self._nbr_idx[lo:hi]

    def bond_satisfaction(self, nd: int) -> tuple[np.ndarray, np.ndarray]:
        """Neighbors of *nd* and each bond's embedded-Ising satisfaction
        ``J_ij (n_i · r)(n_j · r)`` (ref M6) — HALF the energy cost of
        reflecting one endpoint, so the engines' FK activation
        ``1 − e^{−2·sat/T}`` is exactly M6's O(n) bond probability.
        Reflecting an endpoint flips the sign of its ``n·r``, so the
        satisfaction is flip-covariant and signed couplings are valid."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        idx = self._nbr_idx[lo:hi]
        r = self._refl_dir
        proj_nd = float(np.dot(self.s[nd], r))
        return idx, self._nbr_J[lo:hi] * proj_nd * (self.s[idx] @ r)

    def edge_satisfactions(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Per-sweep bond sweep (Swendsen–Wang): ONE fresh reflection
        direction for the whole sweep — every cluster grown and flipped this
        sweep shares it (the embedded-Ising decomposition is per-sweep)."""
        self._draw_reflection()
        proj = self.s @ self._refl_dir
        return (
            self._edge_u,
            self._edge_v,
            self._edge_J * proj[self._edge_u] * proj[self._edge_v],
        )

    def flip_cluster(self, nodes: np.ndarray) -> float:
        """Reflect the cluster's spins about the hyperplane perpendicular to
        the current direction r (``n → n − 2(n·r) r``), in place; return the
        energy change.

        ΔE comes from the boundary bonds only — a joint reflection preserves
        every within-cluster dot product (reflections are orthogonal maps).
        In wolff mode a fresh r is drawn AFTER the flip (one direction per
        cluster, M6); in sw mode the sweep's r is kept (one direction per
        sweep, drawn in :meth:`edge_satisfactions`).
        """
        nodes = np.asarray(nodes, dtype=np.intp)
        r = self._refl_dir
        old = self.s[nodes]
        reflected = old - 2.0 * (old @ r)[:, None] * r[None, :]
        mask = np.zeros(self.N, dtype=np.bool_)
        mask[nodes] = True
        dE = 0.0
        for k, nd in enumerate(nodes):
            lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
            idx = self._nbr_idx[lo:hi]
            out = ~mask[idx]
            J_out = self._nbr_J[lo:hi][out]
            s_out = self.s[idx[out]]
            diff = old[k] - reflected[k]
            dE += float((J_out * (s_out @ diff)).sum())
        self.s[nodes] = reflected
        if self._redraw_dir_on_flip:
            self._draw_reflection()
        return dE

    # --- exchange contract (statsys._thermal.ExchangeSubstrate) ---
    def swap_delta_E(self, u: int, v: int) -> float:
        """ΔE of swapping the spins at *u*, *v*: the two single-move ΔE's
        each count the shared bond as moving to full alignment (−J_uv), but
        a swap leaves the pair's dot product invariant — hence the
        ``+2 J_uv (1 − n_u · n_v)`` correction."""
        dE = self.delta_E(u, self.s[v]) + self.delta_E(v, self.s[u])
        lo, hi = self._nbr_ptr[u], self._nbr_ptr[u + 1]
        row = self._nbr_idx[lo:hi]
        pos = np.flatnonzero(row == v)
        J_uv = float(self._nbr_J[lo:hi][pos[0]]) if len(pos) else 0.0
        return dE + 2.0 * J_uv * (1.0 - float(np.dot(self.s[u], self.s[v])))

    def commit_swap(self, u: int, v: int) -> None:
        su = self.s[u].copy()
        self.s[u] = self.s[v]
        self.s[v] = su

    # ------------------------------------------------------------------
    # Energy (same J as delta_E, by construction)
    # ------------------------------------------------------------------
    def total_energy(self, s: np.ndarray | None = None) -> float:
        """Extensive ``H = −Σ_{(i,j)} J_ij (n_i · n_j)`` (each undirected
        edge once)."""
        self._ensure_couplings()
        state = self.s if s is None else np.asarray(s)
        dots = (state[self._edge_u] * state[self._edge_v]).sum(axis=1)
        return float(-(self._edge_J * dots).sum())

    def energy_intensive(self, s: np.ndarray | None = None) -> float:
        """Per-site energy (the default reporting convention, plan §3b)."""
        return self.total_energy(s) / float(self.N)

    # ------------------------------------------------------------------
    # Lifecycle (RunHostMixin front-door requirements + hooks)
    # ------------------------------------------------------------------
    def init_dynamics(self, custom: Any = None) -> None:
        """Initialise the state and the coupling arrays."""
        self.reset_observables()
        self.init_state(custom)
        self.s = np.asarray(self.s, dtype=np.float64)
        self._ensure_couplings()
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        if getattr(self, "sini", None) is None:
            self.init_dynamics()

    def initialize_run_parameters(
        self, steps: int | None = None, simref: float | None = None
    ) -> None:
        self._set_time_controls(steps=steps, simref=simref)

    def _reset_run_state(self) -> None:
        self._snap_seen = 0
        self._e_running = None

    def _is_py_runlang(self) -> bool:
        return not self.runlang.startswith("C")

    def _record_energy(self) -> None:
        self.observables[HEISENBERG_OBS_ENERGY].append(
            self._e_running / float(self.N)
        )

    def _record_magn(self) -> None:
        self.observables[HEISENBERG_OBS_MAGN].append(self.order_parameter())

    def _record_snapshot(self) -> None:
        snap = self.observables[HEISENBERG_OBS_SNAPSHOTS]
        if snap.streaming:
            if self._snap_seen % self.snapshot_every == 0:
                snap.stream_row(self.s.ravel())
            self._snap_seen += 1
        else:
            snap.append_row(self.s.ravel().copy())

    def _record(self) -> None:
        for record in self._recorders:
            record()

    def _guard_snapshot_size(self) -> None:
        n_rec_est = -(-int(self.steps) // self.snapshot_every)  # ceil
        nbytes = n_rec_est * 3 * self.N * np.dtype(np.float64).itemsize
        if nbytes > HEISENBERG_SNAPSHOT_MAX_BYTES:
            raise ValueError(
                f"snapshot trajectory would be ~{nbytes / 2**20:.0f} MiB "
                f"({n_rec_est} snapshots x 3N={3 * self.N}, float64), "
                f"exceeding HEISENBERG_SNAPSHOT_MAX_BYTES="
                f"{HEISENBERG_SNAPSHOT_MAX_BYTES / 2**20:.0f} MiB. Increase "
                f"snapshot_every (currently {self.snapshot_every}), reduce "
                f"steps, or raise the cap in HeisenbergModel/defaults.py."
            )

    # ------------------------------------------------------------------
    # Vectorized sync backend (np/cu) — VectorSyncEngine hooks
    # ------------------------------------------------------------------
    def _vec_bind(self, xp) -> None:
        """Stage the substrate arrays on the engine's array module
        (zero-copy aliases for np; one host-to-device copy for cu)."""
        self._ensure_couplings()
        if xp is np:
            self._v_s = self.s
            self._v_nbr_J = self._nbr_J
            self._v_nbr_idx = self._nbr_idx
            self._v_nbr_rows = self._nbr_rows
        else:
            self._v_s = xp.asarray(self.s)
            self._v_nbr_J = xp.asarray(self._nbr_J)
            self._v_nbr_idx = xp.asarray(self._nbr_idx)
            self._v_nbr_rows = xp.asarray(self._nbr_rows)

    def _vec_sync_host(self, xp) -> None:
        """Mirror the device state back into ``self.s`` (no-op for np,
        where ``_v_s`` IS ``self.s``)."""
        if xp is not np:
            self.s = xp.asnumpy(self._v_s)

    def _vec_propose(self, cls, xp, rng):
        # ``normalize(n + delta·u)`` with u uniform in [−1, 1]³ — the same
        # kernel (and zero-norm guard: keep the current spin) as
        # propose_local.
        cur = self._v_s[cls]
        perturb = cur + self.delta * rng.uniform(-1.0, 1.0, (cls.shape[0], 3))
        norms = xp.linalg.norm(perturb, axis=1, keepdims=True)
        safe = norms > 0.0
        return xp.where(safe, perturb / xp.where(safe, norms, 1.0), cur)

    def _vec_delta_E(self, cls, proposals, xp):
        """``ΔE_i = (n_i − n'_i) · F_i`` with the molecular field
        ``F_i = Σ_j J_ij n_j`` — the sign matches ``delta_E``."""
        F = xp.zeros((self.N, 3), dtype=xp.float64)
        for c in range(3):
            contrib = self._v_nbr_J * self._v_s[self._v_nbr_idx, c]
            xp.add.at(F[:, c], self._v_nbr_rows, contrib)
        diff = self._v_s[cls] - proposals
        return (diff * F[cls]).sum(axis=1)

    def _vec_commit(self, cls, proposals, accept):
        idx = cls[accept]
        self._v_s[idx] = proposals[accept]

    def _make_vec_engine(self, xp):
        raise NotImplementedError(
            "Scheme classes must implement _make_vec_engine() to enable "
            "the vectorized np/cu backends."
        )

    def _np_check_supported(self) -> None:
        raise NotImplementedError(
            f"{type(self).__name__} has no vectorized (np/cu) backend."
        )

    def _run_np(self, tqdm_on: bool = False, xp=None) -> None:
        """Vectorized sync sampling loop — same record cadence and
        incremental-energy invariant as _sample_py. ``xp`` selects the
        array module (numpy default; cupy for the cu backend)."""
        xp = np if xp is None else xp
        if xp is not np:
            # Device runs draw from xp's global RNG; pin it to the run
            # seed the same way DynSys pins numpy's.
            xp.random.seed(self.seed)
        engine = self._make_vec_engine(xp)
        self._vec_bind(xp)
        sweep = engine.compile_step()
        self._e_running = self.total_energy()
        iterator = (
            tqdm.tqdm(range(self.steps), desc=type(self).__name__)
            if tqdm_on
            else range(self.steps)
        )
        for _ in iterator:
            self._vec_sync_host(xp)
            self._record()
            self._e_running += sweep()
        self._vec_sync_host(xp)

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C: one directory per run)
    # ------------------------------------------------------------------
    def _physics_tokens(self) -> list:
        """Model numerics after p= — overridden per substrate."""
        return [("T", float(self.T))]

    def _axis_tokens(self) -> list:
        """Dynamics axes (rule/move/upd/ord/tf) — leaf hook."""
        return []

    def _name_tokens(self) -> list:
        snap_on = HEISENBERG_OBS_SNAPSHOTS in self._selected_obs
        return [
            ("p", float(self.sg.pflip)),
            *self._physics_tokens(),
            ("ns", int(self.steps)),
            (
                "se",
                int(self.snapshot_every) if snap_on else None,
                HEISENBERG_SNAPSHOT_EVERY_DEFAULT,
            ),
            *self._axis_tokens(),
            ("cn", self.coupling_norm, HEISENBERG_COUPLING_NORM_DEFAULT),
            ("ic", self.ic, "uniform"),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _begin_outputs(self) -> None:
        snap = self.observables[HEISENBERG_OBS_SNAPSHOTS]
        snap.reset_stream()
        self._snap_seen = 0
        if not self.savedisk or not self._is_py_runlang():
            return
        rundir = self._run_output_dir()
        rundir.mkdir(parents=True, exist_ok=True)
        self._write_cfg_sidecar(rundir)
        if HEISENBERG_OBS_ENERGY in self._selected_obs:
            self.observables[HEISENBERG_OBS_ENERGY].set_path(
                rundir / f"{HEISENBERG_ENERGY_FBASE}{BIN}"
            )
        if HEISENBERG_OBS_MAGN in self._selected_obs:
            self.observables[HEISENBERG_OBS_MAGN].set_path(
                rundir / f"{HEISENBERG_MAGN_FBASE}{BIN}"
            )
        if HEISENBERG_OBS_SNAPSHOTS in self._selected_obs:
            snap.set_path(rundir / f"{HEISENBERG_SNAPSHOTS_FBASE}{BIN}")
            self._guard_snapshot_size()
            snap.open_stream()

    def _persist_observables(self) -> None:
        if not self.savedisk or not self._is_py_runlang():
            return
        snap = self.observables[HEISENBERG_OBS_SNAPSHOTS]
        snap.close_stream()
        if (
            HEISENBERG_OBS_SNAPSHOTS in self._selected_obs
            and snap.path is not None
        ):
            snap.finalize_batch(subsample=self.snapshot_every)
        if HEISENBERG_OBS_ENERGY in self._selected_obs:
            self.observables[HEISENBERG_OBS_ENERGY].persist()
        if HEISENBERG_OBS_MAGN in self._selected_obs:
            self.observables[HEISENBERG_OBS_MAGN].persist()

    # ------------------------------------------------------------------
    # Engine plumbing (schemes provide the engine; the base runs it)
    # ------------------------------------------------------------------
    def _make_engine(self):
        raise NotImplementedError(
            "Scheme classes (HeisenbergMetropolis, ...) must implement "
            "_make_engine()."
        )

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """Python-backend sampling loop: record at sweep times 0..steps-1,
        advance one compiled sweep, track the energy incrementally from the
        sweep's ΔE (no per-sweep full recompute)."""
        engine = self._make_engine()
        self._engine = engine  # scheme leaves mirror per-move extras off it
        sweep = engine.compile_step()
        self._e_running = self.total_energy()
        iterator = (
            tqdm.tqdm(range(self.steps), desc=type(self).__name__)
            if tqdm_on
            else range(self.steps)
        )
        for _ in iterator:
            self._record()
            self._e_running += sweep()

    # ------------------------------------------------------------------
    # Native pybind backend hooks (the pb solver polymorphs on these)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        """Setup-time capability gate for ``runlang='pb'`` — the solver calls
        this before any output stream opens. Default: not wired; pb-capable
        scheme leaves override with their own axis validation."""
        raise NotImplementedError(
            f"runlang='pb' is not wired for {type(self).__name__}; use "
            "runlang='py'."
        )

    def _pb_csr(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The symmetric coupling CSR in the dtypes the kernel expects —
        the SAME J as delta_E/total_energy (coupling_norm applied)."""
        self._ensure_couplings()
        return (
            np.ascontiguousarray(self._nbr_idx, dtype=np.int64),
            np.ascontiguousarray(self._nbr_J, dtype=np.float64),
            np.ascontiguousarray(self._nbr_ptr, dtype=np.int64),
        )

    def _run_pb(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """Execute the native kernel (pb-capable leaves override)."""
        raise NotImplementedError(
            f"runlang='pb' is not wired for {type(self).__name__}; use "
            "runlang='py'."
        )

    def make_sweep_fn(self):
        """Frame-producer hook (statsys._frames): one compiled engine sweep."""
        self.check_attribute()
        sweep = self._make_engine().compile_step()

        def _sweep() -> None:
            sweep()

        return _sweep

    def run_py(self, tqdm_on: bool = False, **kw: Any) -> None:
        """Direct Python-backend entry (the registry front door is run())."""
        self.check_attribute()
        self._begin_outputs()
        try:
            self._sample_py(tqdm_on=tqdm_on)
        finally:
            self._persist_observables()

    def run_c(self, **kw: Any) -> None:
        raise NotImplementedError(
            "The C-subprocess backend stays on the legacy HeisenbergModel "
            "class (frozen interface, decision D-B5); use runlang='py'/'pb' "
            "or the legacy HeisenbergModel."
        )
