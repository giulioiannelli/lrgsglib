"""IsingBase — the ±1 Ising substrate on the shared statsys spine.

The Layer-0 model class of the Ising refactor (plan §2/§3): it owns the
Hamiltonian (ONE symmetric coupling matrix ``J`` built once at setup, from
which BOTH ``delta_E`` and the energy observable derive — detailed balance and
ΔE↔E consistency by construction, plan §3b), the external field as a
first-class term, the substrate contract consumed by the shared
``statsys._thermal`` engine, and the ObservableSet lifecycle.

It owns NOTHING scheme-specific: temperature policy, acceptance rule and
update axes live on the thin scheme leaves (``IsingMetropolis``,
``IsingSimulatedAnnealing``, ``IsingParallelTempering``, ``IsingCEM``) —
``IsingCEM`` has no ``T`` at all, which is why ``T`` is not here.

This class is NEW code beside the legacy ``IsingDynamics`` god-object
(strangler pattern, plan §11): the legacy class keeps working untouched; both
dispatch through the same ``"ising"`` solver-registry key, and the solver
polymorphs on the instance.

Hamiltonian (plan §3b, ``k_B = 1``)::

    H = − Σ_{(i,j) ∈ E} J_ij s_i s_j − Σ_i h_i s_i
    ΔE_k(flip) = 2 s_k ( Σ_j J_kj s_j + h_k )

with ``J`` from the signed edge weights via ``coupling_norm``
('raw' | 'sym' | 'avg', defaults.py).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Sequence

import numpy as np
import tqdm

from ...config.const import BIN
from .._c_backend import CBackendMixin
from .._csr import build_coupling_csr
from .._observables import ObservableSet, RowTrajectory, ScalarSeries
from .._run_host import RunHostMixin
from ..BinDynSys import BinDynSys
from .defaults import (
    ISING_COUPLING_NORM_DEFAULT,
    ISING_COUPLING_NORMS,
    ISING_ENERGY_FBASE,
    ISING_MAGN_FBASE,
    ISING_OBS_DEFAULT,
    ISING_OBS_ENERGY,
    ISING_OBS_MAGN,
    ISING_OBS_SNAPSHOTS,
    ISING_SNAPSHOT_EVERY_DEFAULT,
    ISING_SNAPSHOT_MAX_BYTES,
    ISING_SNAPSHOTS_FBASE,
    ISING_SOLVER_NAME,
)

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["IsingBase"]

_KNOWN_OBSERVABLES = (ISING_OBS_ENERGY, ISING_OBS_MAGN, ISING_OBS_SNAPSHOTS)


class IsingBase(RunHostMixin, CBackendMixin, BinDynSys):
    """±1 Ising substrate: symmetric-J Hamiltonian + observables + run host.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings ``w_ij``.
    coupling_norm : {'raw', 'sym', 'avg'}
        How the symmetric ``J`` derives from the edge weights (plan §3b).
    observables : sequence of str
        Which series to record: any of ``'energy'`` (intensive per-site
        energy, tracked incrementally from the sweep ΔE — no per-sweep full
        recompute), ``'magn'`` (per-spin magnetization) and ``'snapshots'``
        ((n_rec, N) int8 trajectory, subsampled by *snapshot_every*).
    snapshot_every : int
        Keep every k-th sweep in the snapshot trajectory.
    **kwargs
        Forwarded to :class:`BinDynSys` (``ic``, ``field``, ``runlang``,
        ``steps``/``simref``, ``seed``, ``savedisk``, ...).
    """

    solver_name = ISING_SOLVER_NAME

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        coupling_norm: str = ISING_COUPLING_NORM_DEFAULT,
        observables: Sequence[str] = ISING_OBS_DEFAULT,
        snapshot_every: int = ISING_SNAPSHOT_EVERY_DEFAULT,
        **kwargs: Any,
    ) -> None:
        if coupling_norm not in ISING_COUPLING_NORMS:
            raise ValueError(
                f"coupling_norm must be one of {ISING_COUPLING_NORMS}; "
                f"got {coupling_norm!r}."
            )
        selected = tuple(observables)
        unknown = [o for o in selected if o not in _KNOWN_OBSERVABLES]
        if unknown:
            raise ValueError(
                f"Unknown observable(s) {unknown}; known: {_KNOWN_OBSERVABLES}."
            )
        # Built BEFORE super().__init__: BinDynSys assigns ``self.magn = []``,
        # which routes through the lazy-view property below and needs the set.
        self.observables = ObservableSet(
            [
                ScalarSeries(ISING_OBS_ENERGY, dtype=np.float64),
                ScalarSeries(ISING_OBS_MAGN, dtype=np.float64),
                RowTrajectory(ISING_OBS_SNAPSHOTS, ncols=sg.N, dtype=np.int8),
            ]
        )
        self._selected_obs = selected
        self.coupling_norm = coupling_norm
        self.snapshot_every = max(1, int(snapshot_every))

        dynpath = getattr(sg, "path_ising", None)
        super().__init__(sg, dynpath=dynpath, **kwargs)

        # Configure-then-run: the recorder list is compiled once here; the
        # per-sweep _record only iterates it (un-selected observables are
        # absent, not skipped by an if).
        recorders = []
        if ISING_OBS_ENERGY in selected:
            recorders.append(self._record_energy)
        if ISING_OBS_MAGN in selected:
            recorders.append(self._record_magn)
        if ISING_OBS_SNAPSHOTS in selected:
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
        return self.observables[ISING_OBS_ENERGY].data

    @ene.setter
    def ene(self, value) -> None:
        self.observables[ISING_OBS_ENERGY].set_values(value)

    @property
    def magn(self) -> list:
        """Per-spin magnetization series m(t) = ⟨s⟩(t)."""
        return self.observables[ISING_OBS_MAGN].data

    @magn.setter
    def magn(self, value) -> None:
        self.observables[ISING_OBS_MAGN].set_values(value)

    @property
    def s_t(self):
        """Spin trajectory: RAM list pre-persist, read-only memmap after."""
        return self.observables[ISING_OBS_SNAPSHOTS].data

    @s_t.setter
    def s_t(self, value) -> None:
        self.observables[ISING_OBS_SNAPSHOTS].set_rows(value)

    # ------------------------------------------------------------------
    # Couplings (ONE symmetric J, built once — plan §3b)
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
        return int(-self.s[nd])

    def local_field(self, nd: int) -> float:
        """``Σ_j J_nd,j s_j + h_nd`` — the field entering ΔE."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        coupling = float(self._nbr_J[lo:hi] @ self.s[self._nbr_idx[lo:hi]])
        return coupling + float(self.field[nd])

    def delta_E(self, nd: int, proposal: int) -> float:
        return (float(self.s[nd]) - float(proposal)) * self.local_field(nd)

    def commit(self, nd: int, proposal: int) -> None:
        self.s[nd] = proposal

    def order_parameter(self) -> float:
        return self.magnetization()

    # --- cluster contract (statsys._thermal.ClusterSubstrate) ---
    def neighbor_indices(self, nd: int) -> np.ndarray:
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        return self._nbr_idx[lo:hi]

    def bond_satisfaction(self, nd: int) -> tuple[np.ndarray, np.ndarray]:
        """Neighbors of *nd* and each bond's satisfaction ``J_ij s_i s_j``
        (> 0 = satisfied = activatable by the FK rule)."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        idx = self._nbr_idx[lo:hi]
        sat = (
            self._nbr_J[lo:hi]
            * float(self.s[nd])
            * self.s[idx].astype(np.float64, copy=False)
        )
        return idx, sat

    def edge_satisfactions(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        s = self.s.astype(np.float64, copy=False)
        sat = self._edge_J * s[self._edge_u] * s[self._edge_v]
        return self._edge_u, self._edge_v, sat

    def flip_cluster(self, nodes: np.ndarray) -> float:
        """Flip *nodes* in place; ΔE comes from the boundary bonds only
        (within-cluster bonds are invariant under a joint flip) plus the
        field term: ``ΔE = 2 Σ_{i∈C} s_i (Σ_{j∉C} J_ij s_j + h_i)``."""
        nodes = np.asarray(nodes, dtype=np.intp)
        mask = np.zeros(self.N, dtype=np.bool_)
        mask[nodes] = True
        s = self.s
        field = np.asarray(self.field, dtype=np.float64)
        dE = 0.0
        for nd in nodes:
            lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
            idx = self._nbr_idx[lo:hi]
            out = ~mask[idx]
            coupling_out = float(self._nbr_J[lo:hi][out] @ s[idx[out]])
            dE += 2.0 * float(s[nd]) * (coupling_out + float(field[nd]))
        self.s[nodes] = -self.s[nodes]
        return dE

    # --- exchange contract (statsys._thermal.ExchangeSubstrate) ---
    def swap_delta_E(self, u: int, v: int) -> float:
        """ΔE of swapping opposite spins at *u*, *v* (= flipping both):
        ``ΔE_u + ΔE_v − 4 J_uv s_u s_v`` — the shared bond is invariant
        under the joint flip, so its two single-flip contributions cancel."""
        dE_u = self.delta_E(u, int(-self.s[u]))
        dE_v = self.delta_E(v, int(-self.s[v]))
        lo, hi = self._nbr_ptr[u], self._nbr_ptr[u + 1]
        row = self._nbr_idx[lo:hi]
        pos = np.flatnonzero(row == v)
        J_uv = float(self._nbr_J[lo:hi][pos[0]]) if len(pos) else 0.0
        return dE_u + dE_v - 4.0 * J_uv * float(self.s[u]) * float(self.s[v])

    def commit_swap(self, u: int, v: int) -> None:
        self.s[u], self.s[v] = self.s[v], self.s[u]

    # ------------------------------------------------------------------
    # Energy (same J as delta_E, by construction)
    # ------------------------------------------------------------------
    def local_fields_all(self) -> np.ndarray:
        """Vectorized ``Σ_j J_kj s_j + h_k`` for every site k."""
        s = self.s.astype(np.float64, copy=False)
        contrib = self._nbr_J * s[self._nbr_idx]
        coupling = np.bincount(
            self._nbr_rows, weights=contrib, minlength=self.N
        )
        return coupling + np.asarray(self.field, dtype=np.float64)

    def all_delta_E(self) -> np.ndarray:
        """ΔE of flipping each site, for the current state."""
        return (
            2.0
            * self.s.astype(np.float64, copy=False)
            * self.local_fields_all()
        )

    def total_energy(self, s: np.ndarray | None = None) -> float:
        """Extensive ``H = −½ Σ_k s_k (J s)_k − h·s`` (CSR is symmetric, so
        the pair sum is half the per-site sum)."""
        self._ensure_couplings()
        state = self.s if s is None else np.asarray(s)
        state = state.astype(np.float64, copy=False)
        contrib = self._nbr_J * state[self._nbr_idx]
        coupling = np.bincount(
            self._nbr_rows, weights=contrib, minlength=self.N
        )
        field = np.asarray(self.field, dtype=np.float64)
        return float(-0.5 * (state @ coupling) - field @ state)

    def energy_intensive(self, s: np.ndarray | None = None) -> float:
        """Per-site energy (the default reporting convention, plan §3b)."""
        return self.total_energy(s) / float(self.N)

    def is_absorbing(self) -> bool:
        """Zero-temperature frozen check on the CURRENT Hamiltonian: no
        downhill flip exists (ties count as mobile unless the scheme's
        ``tie_flip_p`` is 0)."""
        de = self.all_delta_E()
        tie_flip_p = getattr(self, "tie_flip_p", None)
        if tie_flip_p is not None and float(tie_flip_p) <= 0.0:
            return bool(np.all(de >= 0.0))
        return bool(np.all(de > 0.0))

    # ------------------------------------------------------------------
    # Lifecycle (RunHostMixin front-door requirements + hooks)
    # ------------------------------------------------------------------
    def init_dynamics(self, custom: Any = None) -> None:
        """Initialise the state and the coupling arrays (py backend; the C
        export path joins in the native phase)."""
        self.reset_observables()
        self.init_s(custom)
        self.s = np.asarray(self.s, dtype=np.int8)
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
        self.observables[ISING_OBS_ENERGY].append(
            self._e_running / float(self.N)
        )

    def _record_magn(self) -> None:
        self.observables[ISING_OBS_MAGN].append(self.magnetization())

    def _record_snapshot(self) -> None:
        snap = self.observables[ISING_OBS_SNAPSHOTS]
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
        nbytes = n_rec_est * self.N  # int8 => 1 byte/element
        if nbytes > ISING_SNAPSHOT_MAX_BYTES:
            raise ValueError(
                f"snapshot trajectory would be ~{nbytes / 2**20:.0f} MiB "
                f"({n_rec_est} snapshots x N={self.N}), exceeding "
                f"ISING_SNAPSHOT_MAX_BYTES="
                f"{ISING_SNAPSHOT_MAX_BYTES / 2**20:.0f} MiB. Increase "
                f"snapshot_every (currently {self.snapshot_every}), reduce "
                f"steps, or raise the cap in IsingDynamics/defaults.py."
            )

    def _begin_outputs(self) -> None:
        snap = self.observables[ISING_OBS_SNAPSHOTS]
        snap.reset_stream()
        self._snap_seen = 0
        if not self.savedisk or not self._is_py_runlang():
            return
        suf = self.out_suffix
        if ISING_OBS_ENERGY in self._selected_obs:
            self.observables[ISING_OBS_ENERGY].set_path(
                self.dynpath
                / self.sg.get_p_fname(ISING_ENERGY_FBASE, suf, ext=BIN)
            )
        if ISING_OBS_MAGN in self._selected_obs:
            self.observables[ISING_OBS_MAGN].set_path(
                self.dynpath
                / self.sg.get_p_fname(ISING_MAGN_FBASE, suf, ext=BIN)
            )
        if ISING_OBS_SNAPSHOTS in self._selected_obs:
            snap.set_path(
                self.dynpath
                / self.sg.get_p_fname(ISING_SNAPSHOTS_FBASE, suf, ext=BIN)
            )
            self._guard_snapshot_size()
            snap.open_stream()

    def _persist_observables(self) -> None:
        if not self.savedisk or not self._is_py_runlang():
            return
        snap = self.observables[ISING_OBS_SNAPSHOTS]
        snap.close_stream()
        if ISING_OBS_SNAPSHOTS in self._selected_obs and snap.path is not None:
            snap.finalize_batch(subsample=self.snapshot_every)
        if ISING_OBS_ENERGY in self._selected_obs:
            self.observables[ISING_OBS_ENERGY].persist()
        if ISING_OBS_MAGN in self._selected_obs:
            self.observables[ISING_OBS_MAGN].persist()

    # ------------------------------------------------------------------
    # Engine plumbing (schemes provide the engine; the base runs it)
    # ------------------------------------------------------------------
    def _make_engine(self):
        raise NotImplementedError(
            "Scheme classes (IsingMetropolis, ...) must implement _make_engine()."
        )

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """Python-backend sampling loop: record at sweep times 0..steps-1
        (the Voter/CP cadence), advance one compiled sweep, track the energy
        incrementally from the sweep's ΔE (no per-sweep full recompute)."""
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

    def _pb_check_kinetics(
        self,
        *,
        rule: str,
        order: str,
        tie_flip_p: float | None,
        zero_T: bool = False,
    ) -> None:
        """Shared guards for the native single-site kernel
        (``LRGSG_rbim.c::glauberMetropolis``): whatever the compiled kernel
        cannot represent is a hard capability error, never a silent
        approximation (plan invariant #3).

        The kernel implements the Metropolis acceptance only, visits sites
        in 'sequential' or 'asynchronous' order, bakes its tie policy
        (``e^0 = 1`` at T > 0; frozen at T = 0), and — crucially — records
        the pairwise-only energy while using the field in ΔE (the legacy
        convention), so a nonzero field would silently corrupt the §3b
        energy trace.
        """
        if rule != "metropolis":
            raise NotImplementedError(
                "runlang='pb': the native kernel implements the metropolis "
                f"acceptance only; rule={rule!r} is python-only."
            )
        if order == "permutation":
            raise NotImplementedError(
                "runlang='pb': the native kernel visits sites in 'random' "
                "or 'typewriter' order; order='permutation' is python-only."
            )
        expected_tie = 0.0 if zero_T else 1.0
        if tie_flip_p is not None and float(tie_flip_p) != expected_tie:
            raise NotImplementedError(
                "runlang='pb': the native kernel bakes its ΔE=0 policy "
                f"(accept with p={expected_tie} at this temperature); "
                f"tie_flip_p={tie_flip_p} is python-only."
            )
        if np.any(np.asarray(self.field, dtype=np.float64) != 0.0):
            raise NotImplementedError(
                "runlang='pb': the native kernel records the pairwise-only "
                "energy (field enters ΔE but not the recorded trace — the "
                "legacy convention), which would break the field-first-class "
                "energy invariant; use h = 0 or runlang='py'."
            )
        if ISING_OBS_SNAPSHOTS in self._selected_obs:
            raise NotImplementedError(
                "runlang='pb': the native kernel returns energy/magn traces "
                "and the final state only; deselect 'snapshots' or use "
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
            "The C backend for the new Ising scheme classes lands in Phase 2 "
            "(native backends); use runlang='py' or the legacy IsingDynamics."
        )
