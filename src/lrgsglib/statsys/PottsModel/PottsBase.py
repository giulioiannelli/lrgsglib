"""PottsBase — the q-state Potts substrate on the shared statsys spine.

The Layer-0 model class of the Potts modernization (statsys unification
D-B8/D-B9; the second consumer of the Ising-refactor architecture, proving
the write-once claim of the shared engine): it owns the Hamiltonian (ONE
symmetric coupling matrix ``J`` built once by the shared
``statsys._csr.build_coupling_csr``, from which BOTH ``delta_E`` and the
energy observable derive — ΔE↔E consistency by construction, plan §3b), the
substrate contracts consumed by the shared ``statsys._thermal`` engines
(single-site + cluster + exchange), and the ObservableSet lifecycle.

It owns NOTHING scheme-specific — temperature policy, acceptance rule and
update axes live on the thin scheme leaves (``PottsMetropolis``, ...).

This class is NEW code beside the legacy ``PottsModel`` (strangler pattern):
the legacy class keeps working untouched; both dispatch through the same
``"potts"`` solver-registry key, and the solver polymorphs on the instance.

Hamiltonian (ref M10, Wu 1982; ``k_B = 1``)::

    H = − Σ_{(i,j) ∈ E} J_ij δ(σ_i, σ_j),      σ_i ∈ {0, ..., q−1}
    ΔE(i: a → b) = Σ_j J_ij [δ(a, σ_j) − δ(b, σ_j)]

with ``J`` from the signed edge weights via ``coupling_norm``. The Potts
ordering-field term is NOT wired: the external field must be zero (hard
capability error, invariant #3). For ``q = 2`` under ``σ → s = 2σ − 1`` this
is the Ising model at halved couplings (``δ = (1 + s_i s_j)/2``).
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
from ..VecDynSys import VecDynSys
from .defaults import (
    POTTS_COUPLING_NORM_DEFAULT,
    POTTS_COUPLING_NORMS,
    POTTS_DYN_SUBDIR,
    POTTS_ENERGY_FBASE,
    POTTS_MAGN_FBASE,
    POTTS_OBS_DEFAULT,
    POTTS_OBS_ENERGY,
    POTTS_OBS_MAGN,
    POTTS_OBS_SNAPSHOTS,
    POTTS_Q_DEFAULT,
    POTTS_SNAPSHOT_EVERY_DEFAULT,
    POTTS_SNAPSHOT_MAX_BYTES,
    POTTS_SNAPSHOTS_FBASE,
    POTTS_SOLVER_NAME,
    POTTS_T_DEFAULT,
)

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["PottsBase"]

_KNOWN_OBSERVABLES = (POTTS_OBS_ENERGY, POTTS_OBS_MAGN, POTTS_OBS_SNAPSHOTS)


class PottsBase(RunHostMixin, CBackendMixin, VecDynSys):
    """q-state Potts substrate: symmetric-J Hamiltonian + observables + run host.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings ``w_ij``.
    q : int
        Number of Potts states (>= 2).
    T : float
        Temperature (``k_B = 1``; stored by ``VecDynSys``, the policy lives
        on the scheme leaves).
    coupling_norm : {'raw', 'sym', 'avg'}
        How the symmetric ``J`` derives from the edge weights (shared
        vocabulary, ``statsys._csr``).
    observables : sequence of str
        Which series to record: any of ``'energy'`` (intensive per-site
        energy, tracked incrementally from the sweep ΔE), ``'magn'`` (the
        Potts order parameter ``(q·max_frac − 1)/(q − 1)``, ref M10 — the
        name keeps the cross-model series convention) and ``'snapshots'``
        ((n_rec, N) int32 label trajectory, subsampled by *snapshot_every*).
    snapshot_every : int
        Keep every k-th sweep in the snapshot trajectory.
    **kwargs
        Forwarded to :class:`VecDynSys` (``ic``, ``runlang``,
        ``steps``/``simref``, ``seed``, ``savedisk``, ...).
    """

    solver_name = POTTS_SOLVER_NAME

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        q: int = POTTS_Q_DEFAULT,
        T: float = POTTS_T_DEFAULT,
        coupling_norm: str = POTTS_COUPLING_NORM_DEFAULT,
        observables: Sequence[str] = POTTS_OBS_DEFAULT,
        snapshot_every: int = POTTS_SNAPSHOT_EVERY_DEFAULT,
        **kwargs: Any,
    ) -> None:
        if int(q) < 2:
            raise ValueError(f"q must be >= 2; got {q}.")
        if coupling_norm not in POTTS_COUPLING_NORMS:
            raise ValueError(
                f"coupling_norm must be one of {POTTS_COUPLING_NORMS}; "
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
                ScalarSeries(POTTS_OBS_ENERGY, dtype=np.float64),
                ScalarSeries(POTTS_OBS_MAGN, dtype=np.float64),
                RowTrajectory(POTTS_OBS_SNAPSHOTS, ncols=sg.N, dtype=np.int32),
            ]
        )
        self._selected_obs = selected
        self.coupling_norm = coupling_norm
        self.snapshot_every = max(1, int(snapshot_every))

        dynpath = getattr(sg, "path_data", None)
        if dynpath is not None:
            dynpath = Path(dynpath) / POTTS_DYN_SUBDIR
        super().__init__(
            sg,
            q=int(q),
            discrete=True,
            T=float(T),
            dynpath=dynpath,
            **kwargs,
        )
        # The Potts ordering-field term (a per-(site,state) coupling) is not
        # wired; DynSys's per-site field vector has no Potts meaning yet.
        if np.any(np.asarray(self.field, dtype=np.float64) != 0.0):
            raise NotImplementedError(
                "The Potts external-field (ordering-field) term is not "
                "wired; the field must be zero."
            )

        # Configure-then-run: the recorder list is compiled once here; the
        # per-sweep _record only iterates it.
        recorders = []
        if POTTS_OBS_ENERGY in selected:
            recorders.append(self._record_energy)
        if POTTS_OBS_MAGN in selected:
            recorders.append(self._record_magn)
        if POTTS_OBS_SNAPSHOTS in selected:
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
        # Running extensive energy, advanced by the engine's per-sweep ΔE.
        self._e_running: float | None = None
        # The engine of the last run (per-move extras live on it).
        self._engine: Any = None
        self.sini: np.ndarray | None = None
        self.reset_observables()

    # ------------------------------------------------------------------
    # Disk-backed observables (lazy views over the ObservableSet)
    # ------------------------------------------------------------------
    @property
    def ene(self) -> list:
        """Intensive per-site energy series e(t) = E(t)/N."""
        return self.observables[POTTS_OBS_ENERGY].data

    @ene.setter
    def ene(self, value) -> None:
        self.observables[POTTS_OBS_ENERGY].set_values(value)

    @property
    def magn(self) -> list:
        """Order-parameter series m(t) = (q·max_frac − 1)/(q − 1)."""
        return self.observables[POTTS_OBS_MAGN].data

    @magn.setter
    def magn(self, value) -> None:
        self.observables[POTTS_OBS_MAGN].set_values(value)

    @property
    def s_t(self):
        """Label trajectory: RAM list pre-persist, read-only memmap after."""
        return self.observables[POTTS_OBS_SNAPSHOTS].data

    @s_t.setter
    def s_t(self, value) -> None:
        self.observables[POTTS_OBS_SNAPSHOTS].set_rows(value)

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
        # Each undirected edge once (u < v), for the SW bond sweep.
        self._edge_u = c.edge_u
        self._edge_v = c.edge_v
        self._edge_J = c.edge_J

    # ------------------------------------------------------------------
    # Substrate contract (consumed by statsys._thermal.ThermalEngine)
    # ------------------------------------------------------------------
    def propose_local(self, nd: int) -> int:
        """A uniform random label different from the current one (for q = 2
        this is deterministic — the only other label)."""
        r = int(np.random.randint(self.q - 1))
        return r + 1 if r >= int(self.s[nd]) else r

    def _label_field(self, nd: int, label: int) -> float:
        """``Σ_j J_nd,j δ(label, σ_j)`` — the bond alignment site *nd* would
        have with *label* (its energy contribution is minus this)."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        nbr_s = self.s[self._nbr_idx[lo:hi]]
        return float(self._nbr_J[lo:hi][nbr_s == label].sum())

    def delta_E(self, nd: int, proposal: int) -> float:
        return self._label_field(nd, int(self.s[nd])) - self._label_field(
            nd, int(proposal)
        )

    def commit(self, nd: int, proposal: int) -> None:
        self.s[nd] = proposal

    def order_parameter(self, s: np.ndarray | None = None) -> float:
        """Potts order parameter ``(q·max_frac − 1)/(q − 1)`` (ref M10):
        1 in a monochromatic state, 0 for equipopulated labels."""
        state = self.s if s is None else np.asarray(s)
        counts = np.bincount(state.astype(np.int64), minlength=self.q)
        max_frac = counts.max() / float(self.N)
        return (self.q * max_frac - 1.0) / (self.q - 1.0)

    def per_state_density(self, s: np.ndarray | None = None) -> np.ndarray:
        """Fraction of sites carrying each label."""
        state = self.s if s is None else np.asarray(s)
        counts = np.bincount(state.astype(np.int64), minlength=self.q)
        return counts / float(self.N)

    # --- cluster contract (statsys._thermal.ClusterSubstrate) ---
    def neighbor_indices(self, nd: int) -> np.ndarray:
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        return self._nbr_idx[lo:hi]

    def bond_satisfaction(self, nd: int) -> tuple[np.ndarray, np.ndarray]:
        """Neighbors of *nd* and each bond's satisfaction
        ``J_ij δ(σ_i, σ_j)/2`` — HALF the cost of breaking a same-label bond,
        so the engines' FK activation ``1 − e^{−2·sat/T}`` is ref M5's
        ``1 − e^{−βJ}``."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        idx = self._nbr_idx[lo:hi]
        same = (self.s[idx] == self.s[nd]).astype(np.float64)
        return idx, 0.5 * self._nbr_J[lo:hi] * same

    def edge_satisfactions(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        same = (self.s[self._edge_u] == self.s[self._edge_v]).astype(np.float64)
        return self._edge_u, self._edge_v, 0.5 * self._edge_J * same

    def flip_cluster(self, nodes: np.ndarray) -> float:
        """Relabel the (same-label) cluster *nodes* to a uniform random OTHER
        label, in place; return the energy change.

        ΔE comes from the boundary bonds only — within-cluster bonds join two
        same-label sites before AND after the joint relabel, so their δ is
        invariant: ``ΔE = Σ_{i∈C, j∉C} J_ij [δ(σ, σ_j) − δ(τ, σ_j)]``. Valid
        for non-negative couplings only (the scheme guards): an
        antiferromagnetic bond's satisfied state is an UNLIKE pair, which a
        joint relabel does not preserve, so the Ising flip-covariance
        argument does not carry over.
        """
        nodes = np.asarray(nodes, dtype=np.intp)
        sigma = int(self.s[nodes[0]])
        r = int(np.random.randint(self.q - 1))
        tau = r + 1 if r >= sigma else r
        mask = np.zeros(self.N, dtype=np.bool_)
        mask[nodes] = True
        dE = 0.0
        for nd in nodes:
            lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
            idx = self._nbr_idx[lo:hi]
            out = ~mask[idx]
            J_out = self._nbr_J[lo:hi][out]
            s_out = self.s[idx[out]]
            dE += float(J_out[s_out == sigma].sum()) - float(
                J_out[s_out == tau].sum()
            )
        self.s[nodes] = tau
        return dE

    # --- exchange contract (statsys._thermal.ExchangeSubstrate) ---
    def swap_delta_E(self, u: int, v: int) -> float:
        """ΔE of swapping the (different) labels at *u*, *v*: the two
        single-move ΔE's each count the shared bond as broken (−J_uv), but a
        swap leaves the pair unlike — the bond term is invariant — hence the
        ``+2 J_uv`` correction."""
        dE = self.delta_E(u, int(self.s[v])) + self.delta_E(v, int(self.s[u]))
        lo, hi = self._nbr_ptr[u], self._nbr_ptr[u + 1]
        row = self._nbr_idx[lo:hi]
        pos = np.flatnonzero(row == v)
        J_uv = float(self._nbr_J[lo:hi][pos[0]]) if len(pos) else 0.0
        return dE + 2.0 * J_uv

    def commit_swap(self, u: int, v: int) -> None:
        self.s[u], self.s[v] = self.s[v], self.s[u]

    # ------------------------------------------------------------------
    # Energy (same J as delta_E, by construction)
    # ------------------------------------------------------------------
    def label_fields_all(self, s: np.ndarray | None = None) -> np.ndarray:
        """Vectorized ``F[i, a] = Σ_j J_ij δ(a, σ_j)`` for every site and
        label — site *i*'s bond alignment with each candidate label."""
        self._ensure_couplings()
        state = self.s if s is None else np.asarray(s)
        F = np.zeros((self.N, self.q), dtype=np.float64)
        np.add.at(
            F,
            (self._nbr_rows, state[self._nbr_idx].astype(np.int64)),
            self._nbr_J,
        )
        return F

    def all_delta_E(self) -> np.ndarray:
        """(N, q) matrix of ΔE for moving each site to each label (0 on the
        current label), for the current state."""
        F = self.label_fields_all()
        own = F[np.arange(self.N), self.s.astype(np.int64)]
        return own[:, None] - F

    def total_energy(self, s: np.ndarray | None = None) -> float:
        """Extensive ``H = −Σ_{(i,j)} J_ij δ(σ_i, σ_j)`` (the CSR is
        symmetric, so the pair sum is half the per-site alignment sum)."""
        state = self.s if s is None else np.asarray(s)
        F = self.label_fields_all(state)
        own = F[np.arange(self.N), state.astype(np.int64)]
        return float(-0.5 * own.sum())

    def energy_intensive(self, s: np.ndarray | None = None) -> float:
        """Per-site energy (the default reporting convention, plan §3b)."""
        return self.total_energy(s) / float(self.N)

    def is_absorbing(self) -> bool:
        """Zero-temperature frozen check on the CURRENT Hamiltonian: no
        downhill label change exists (ties count as mobile unless the
        scheme's ``tie_flip_p`` is 0)."""
        de = self.all_delta_E()
        de[np.arange(self.N), self.s.astype(np.int64)] = np.inf  # self-moves
        tie_flip_p = getattr(self, "tie_flip_p", None)
        if tie_flip_p is not None and float(tie_flip_p) <= 0.0:
            return bool(np.all(de >= 0.0))
        return bool(np.all(de > 0.0))

    # ------------------------------------------------------------------
    # Lifecycle (RunHostMixin front-door requirements + hooks)
    # ------------------------------------------------------------------
    def init_dynamics(self, custom: Any = None) -> None:
        """Initialise the state and the coupling arrays."""
        self.reset_observables()
        self.init_state(custom)
        self.s = np.asarray(self.s, dtype=np.int32)
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
        self.observables[POTTS_OBS_ENERGY].append(
            self._e_running / float(self.N)
        )

    def _record_magn(self) -> None:
        self.observables[POTTS_OBS_MAGN].append(self.order_parameter())

    def _record_snapshot(self) -> None:
        snap = self.observables[POTTS_OBS_SNAPSHOTS]
        if snap.streaming:
            if self._snap_seen % self.snapshot_every == 0:
                snap.stream_row(self.s)
            self._snap_seen += 1
        else:
            snap.append_row(self.s.copy())

    def _record(self) -> None:
        for record in self._recorders:
            record()

    def _guard_snapshot_size(self) -> None:
        n_rec_est = -(-int(self.steps) // self.snapshot_every)  # ceil
        nbytes = n_rec_est * self.N * np.dtype(np.int32).itemsize
        if nbytes > POTTS_SNAPSHOT_MAX_BYTES:
            raise ValueError(
                f"snapshot trajectory would be ~{nbytes / 2**20:.0f} MiB "
                f"({n_rec_est} snapshots x N={self.N}, int32), exceeding "
                f"POTTS_SNAPSHOT_MAX_BYTES="
                f"{POTTS_SNAPSHOT_MAX_BYTES / 2**20:.0f} MiB. Increase "
                f"snapshot_every (currently {self.snapshot_every}), reduce "
                f"steps, or raise the cap in PottsModel/defaults.py."
            )

    def _begin_outputs(self) -> None:
        snap = self.observables[POTTS_OBS_SNAPSHOTS]
        snap.reset_stream()
        self._snap_seen = 0
        if not self.savedisk or not self._is_py_runlang():
            return
        suf = self.out_suffix
        if POTTS_OBS_ENERGY in self._selected_obs:
            self.observables[POTTS_OBS_ENERGY].set_path(
                self.dynpath
                / self.sg.get_p_fname(POTTS_ENERGY_FBASE, suf, ext=BIN)
            )
        if POTTS_OBS_MAGN in self._selected_obs:
            self.observables[POTTS_OBS_MAGN].set_path(
                self.dynpath
                / self.sg.get_p_fname(POTTS_MAGN_FBASE, suf, ext=BIN)
            )
        if POTTS_OBS_SNAPSHOTS in self._selected_obs:
            snap.set_path(
                self.dynpath
                / self.sg.get_p_fname(POTTS_SNAPSHOTS_FBASE, suf, ext=BIN)
            )
            self._guard_snapshot_size()
            snap.open_stream()

    def _persist_observables(self) -> None:
        if not self.savedisk or not self._is_py_runlang():
            return
        snap = self.observables[POTTS_OBS_SNAPSHOTS]
        snap.close_stream()
        if POTTS_OBS_SNAPSHOTS in self._selected_obs and snap.path is not None:
            snap.finalize_batch(subsample=self.snapshot_every)
        if POTTS_OBS_ENERGY in self._selected_obs:
            self.observables[POTTS_OBS_ENERGY].persist()
        if POTTS_OBS_MAGN in self._selected_obs:
            self.observables[POTTS_OBS_MAGN].persist()

    # ------------------------------------------------------------------
    # Engine plumbing (schemes provide the engine; the base runs it)
    # ------------------------------------------------------------------
    def _make_engine(self):
        raise NotImplementedError(
            "Scheme classes (PottsMetropolis, ...) must implement _make_engine()."
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
            "The C-subprocess backend stays on the legacy PottsModel class "
            "(frozen interface, decision D-B5); use runlang='py'/'pb' or the "
            "legacy PottsModel."
        )
