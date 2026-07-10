"""MultiSpeciesBase — k coupled Potts components on the shared statsys spine.

The Layer-0 model class of the multi-species modernization (statsys
unification D-B8/D-B9; the last member of the spin quad): it owns the
Hamiltonian (ONE symmetric coupling matrix ``J`` built once by the shared
``statsys._csr.build_coupling_csr``, from which BOTH ``delta_E`` and the
energy observable derive — ΔE↔E consistency by construction, plan §3b), the
substrate contracts consumed by the shared ``statsys._thermal`` engines
(single-site + exchange), and the ObservableSet lifecycle.

It owns NOTHING scheme-specific — temperature policy, acceptance rule and
update axes live on the thin scheme leaves (``MultiSpeciesMetropolis``, ...).

This class is NEW code beside the legacy ``MultiSpeciesModel`` (strangler
pattern): the legacy class keeps working untouched; both dispatch through the
same ``"multispec"`` solver-registry key, and the solver polymorphs on the
instance.

Hamiltonian (k coupled Potts components; ``k_B = 1``)::

    H = − Σ_{(i,j) ∈ E} Σ_{s1,s2} J_ij M[s1,s2] δ(σ_i^{s1}, σ_j^{s2})

with ``J`` from the signed edge weights via ``coupling_norm`` and ``M`` the
(k, k) species-interaction matrix. ``M`` must be SYMMETRIC (hard capability
error): in the pair term site *i* always occupies the first species slot, so
only a symmetric ``M`` makes the single-component ΔE independent of which
slot the updated site sits in — i.e. consistent with the Hamiltonian.
``k = 1`` with ``M = 1`` is exactly the q-state Potts model (ref M10).
The external field is NOT wired (must be zero, invariant #3).
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
    MULTISPEC_COUPLING_NORM_DEFAULT,
    MULTISPEC_COUPLING_NORMS,
    MULTISPEC_DYN_SUBDIR,
    MULTISPEC_ENERGY_FBASE,
    MULTISPEC_MAGN_FBASE,
    MULTISPEC_OBS_DEFAULT,
    MULTISPEC_OBS_ENERGY,
    MULTISPEC_OBS_MAGN,
    MULTISPEC_OBS_SNAPSHOTS,
    MULTISPEC_Q_DEFAULT,
    MULTISPEC_SNAPSHOT_EVERY_DEFAULT,
    MULTISPEC_SNAPSHOT_MAX_BYTES,
    MULTISPEC_SNAPSHOTS_FBASE,
    MULTISPEC_SOLVER_NAME,
    MULTISPEC_SPECIES_DEFAULT,
    MULTISPEC_T_DEFAULT,
)

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["MultiSpeciesBase"]

_KNOWN_OBSERVABLES = (
    MULTISPEC_OBS_ENERGY,
    MULTISPEC_OBS_MAGN,
    MULTISPEC_OBS_SNAPSHOTS,
)


class MultiSpeciesBase(RunHostMixin, CBackendMixin, VecDynSys):
    """k-species substrate: symmetric-J Hamiltonian + observables + run host.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings ``w_ij``.
    species : int
        Number of coupled Potts components per site (k >= 1).
    q_per_species : int or sequence of int
        Number of states per component (each >= 2); a scalar applies to
        every component.
    T : float
        Temperature (``k_B = 1``; stored by ``VecDynSys``, the policy lives
        on the scheme leaves).
    interaction_matrix : ndarray, optional
        The (k, k) species-interaction matrix ``M`` — SYMMETRIC only (see
        the module docstring). Defaults to the identity (each species is an
        independent Potts layer on the same graph).
    coupling_norm : {'raw', 'sym', 'avg'}
        How the symmetric ``J`` derives from the edge weights (shared
        vocabulary, ``statsys._csr``).
    observables : sequence of str
        Which series to record: any of ``'energy'`` (intensive per-site
        energy, tracked incrementally from the sweep ΔE), ``'magn'`` (the
        species-AVERAGED Potts order parameter
        ``mean_k (q_k·max_frac_k − 1)/(q_k − 1)``, ref M10) and
        ``'snapshots'`` ((n_rec, kN) int32 trajectory — each row is the
        flattened (N, k) state, subsampled by *snapshot_every*).
    snapshot_every : int
        Keep every k-th sweep in the snapshot trajectory.
    **kwargs
        Forwarded to :class:`VecDynSys` (``ic``, ``runlang``,
        ``steps``/``simref``, ``seed``, ``savedisk``, ...).
    """

    solver_name = MULTISPEC_SOLVER_NAME

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        species: int = MULTISPEC_SPECIES_DEFAULT,
        q_per_species: int | Sequence[int] = MULTISPEC_Q_DEFAULT,
        T: float = MULTISPEC_T_DEFAULT,
        interaction_matrix: np.ndarray | None = None,
        coupling_norm: str = MULTISPEC_COUPLING_NORM_DEFAULT,
        observables: Sequence[str] = MULTISPEC_OBS_DEFAULT,
        snapshot_every: int = MULTISPEC_SNAPSHOT_EVERY_DEFAULT,
        **kwargs: Any,
    ) -> None:
        if int(species) < 1:
            raise ValueError(f"species must be >= 1; got {species}.")
        self.species = int(species)
        if isinstance(q_per_species, (int, np.integer)):
            self._q_per_species = [int(q_per_species)] * self.species
        else:
            self._q_per_species = [int(q) for q in q_per_species]
            if len(self._q_per_species) != self.species:
                raise ValueError(
                    "q_per_species length must match the species count; "
                    f"got {len(self._q_per_species)} != {self.species}."
                )
        if any(q < 2 for q in self._q_per_species):
            raise ValueError(
                f"every q_per_species must be >= 2; got {self._q_per_species}."
            )
        if interaction_matrix is None:
            M = np.eye(self.species, dtype=np.float64)
        else:
            M = np.asarray(interaction_matrix, dtype=np.float64)
            if M.shape != (self.species, self.species):
                raise ValueError(
                    f"interaction_matrix shape {M.shape} != "
                    f"({self.species}, {self.species})."
                )
            if not np.allclose(M, M.T):
                raise NotImplementedError(
                    "Only a SYMMETRIC species-interaction matrix is wired: "
                    "in the pair term the lower-index site always occupies "
                    "the first species slot, so an asymmetric M makes the "
                    "single-component ΔE inconsistent with the Hamiltonian."
                )
        self.interaction_matrix = M
        if coupling_norm not in MULTISPEC_COUPLING_NORMS:
            raise ValueError(
                f"coupling_norm must be one of {MULTISPEC_COUPLING_NORMS}; "
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
                ScalarSeries(MULTISPEC_OBS_ENERGY, dtype=np.float64),
                ScalarSeries(MULTISPEC_OBS_MAGN, dtype=np.float64),
                # One row per record = the flattened (N, k) state.
                RowTrajectory(
                    MULTISPEC_OBS_SNAPSHOTS,
                    ncols=self.species * sg.N,
                    dtype=np.int32,
                ),
            ]
        )
        self._selected_obs = selected
        self.coupling_norm = coupling_norm
        self.snapshot_every = max(1, int(snapshot_every))

        # Graph-anchored dynamics tree (data/<graph>/<model>/N=...),
        # like path_ising; legacy classes keep the old flat location.
        dynpath = getattr(sg, "path_multi_species", None)
        if dynpath is None:
            base = getattr(sg, "path_data", None)
            if base is not None:
                dynpath = Path(base) / MULTISPEC_DYN_SUBDIR
        super().__init__(
            sg,
            q=max(self._q_per_species),
            discrete=True,
            T=float(T),
            dynpath=dynpath,
            **kwargs,
        )
        # A per-(site, component, label) ordering field is not wired;
        # DynSys's per-site field vector has no multi-species meaning yet.
        if np.any(np.asarray(self.field, dtype=np.float64) != 0.0):
            raise NotImplementedError(
                "The multi-species external-field term is not wired; the "
                "field must be zero."
            )

        # Configure-then-run: the recorder list is compiled once here; the
        # per-sweep _record only iterates it.
        recorders = []
        if MULTISPEC_OBS_ENERGY in selected:
            recorders.append(self._record_energy)
        if MULTISPEC_OBS_MAGN in selected:
            recorders.append(self._record_magn)
        if MULTISPEC_OBS_SNAPSHOTS in selected:
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

    @property
    def _state_shape(self) -> tuple[int, ...]:
        return (self.N, self.species)

    # ------------------------------------------------------------------
    # State initialisation (one label per component)
    # ------------------------------------------------------------------
    def init_state(self, custom: Any = None) -> None:
        """Initialise labels: 'uniform' draws each component uniformly in
        its own [0, q_k); 'zero' sets everything to label 0; 'custom' takes
        an (N, k) integer array (labels validated per component)."""
        match self.ic:
            case "uniform" | "random" | "rand":
                self.s = np.zeros((self.N, self.species), dtype=np.int32)
                for k in range(self.species):
                    self.s[:, k] = np.random.randint(
                        0, self._q_per_species[k], size=self.N
                    )
            case "zero" | "zeros":
                self.s = np.zeros((self.N, self.species), dtype=np.int32)
            case "custom":
                if custom is None:
                    raise ValueError("Provide a custom state array.")
                self.s = np.asarray(custom, dtype=np.int32).copy()
                if self.s.shape != (self.N, self.species):
                    raise ValueError(
                        f"Shape mismatch: {self.s.shape} != "
                        f"({self.N}, {self.species})"
                    )
                for k in range(self.species):
                    q_k = self._q_per_species[k]
                    if ((self.s[:, k] < 0) | (self.s[:, k] >= q_k)).any():
                        raise ValueError(
                            f"custom labels of species {k} must lie in "
                            f"[0, {q_k})."
                        )
            case _:
                raise ValueError(f"Unknown IC: {self.ic!r}")

    # ------------------------------------------------------------------
    # Disk-backed observables (lazy views over the ObservableSet)
    # ------------------------------------------------------------------
    @property
    def ene(self) -> list:
        """Intensive per-site energy series e(t) = E(t)/N."""
        return self.observables[MULTISPEC_OBS_ENERGY].data

    @ene.setter
    def ene(self, value) -> None:
        self.observables[MULTISPEC_OBS_ENERGY].set_values(value)

    @property
    def magn(self) -> list:
        """Species-averaged order-parameter series."""
        return self.observables[MULTISPEC_OBS_MAGN].data

    @magn.setter
    def magn(self, value) -> None:
        self.observables[MULTISPEC_OBS_MAGN].set_values(value)

    @property
    def s_t(self):
        """Label trajectory (rows = flattened (N, k) states): RAM list
        pre-persist, read-only memmap after."""
        return self.observables[MULTISPEC_OBS_SNAPSHOTS].data

    @s_t.setter
    def s_t(self, value) -> None:
        self.observables[MULTISPEC_OBS_SNAPSHOTS].set_rows(value)

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
        # Each undirected edge once (u < v), for the energy sum.
        self._edge_u = c.edge_u
        self._edge_v = c.edge_v
        self._edge_J = c.edge_J

    # ------------------------------------------------------------------
    # Substrate contract (consumed by statsys._thermal.ThermalEngine)
    # ------------------------------------------------------------------
    def propose_local(self, nd: int) -> tuple[int, int]:
        """A uniform random component and a uniform random OTHER label for
        it — the historical (component, label) proposal, symmetric."""
        comp = int(np.random.randint(self.species))
        q_k = self._q_per_species[comp]
        r = int(np.random.randint(q_k - 1))
        cur = int(self.s[nd, comp])
        return comp, (r + 1 if r >= cur else r)

    def _component_field(self, nd: int, comp: int, label: int) -> float:
        """``Σ_j J_nd,j Σ_{s2} M[comp,s2] δ(label, σ_j^{s2})`` — the bond
        alignment component *comp* of site *nd* would have with *label*."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        J = self._nbr_J[lo:hi]
        nbr_s = self.s[self._nbr_idx[lo:hi]]
        # per-species alignment weights, then the interaction row folds them
        per_species = (J[:, None] * (nbr_s == label)).sum(axis=0)
        return float(self.interaction_matrix[comp] @ per_species)

    def delta_E(self, nd: int, proposal: tuple[int, int]) -> float:
        comp, label = proposal
        cur = int(self.s[nd, comp])
        return self._component_field(nd, comp, cur) - self._component_field(
            nd, comp, label
        )

    def commit(self, nd: int, proposal: tuple[int, int]) -> None:
        comp, label = proposal
        self.s[nd, comp] = label

    def order_parameter(self, s: np.ndarray | None = None) -> float:
        """Species-averaged Potts order parameter
        ``mean_k (q_k·max_frac_k − 1)/(q_k − 1)`` (ref M10): 1 when every
        species is monochromatic, 0 for equipopulated labels."""
        state = self.s if s is None else np.asarray(s)
        vals = []
        for k in range(self.species):
            q_k = self._q_per_species[k]
            counts = np.bincount(state[:, k].astype(np.int64), minlength=q_k)
            max_frac = counts.max() / float(self.N)
            vals.append((q_k * max_frac - 1.0) / (q_k - 1.0))
        return float(np.mean(vals))

    def per_species_density(
        self, s: np.ndarray | None = None
    ) -> list[np.ndarray]:
        """Fraction of sites carrying each label, per species."""
        state = self.s if s is None else np.asarray(s)
        return [
            np.bincount(
                state[:, k].astype(np.int64),
                minlength=self._q_per_species[k],
            )
            / float(self.N)
            for k in range(self.species)
        ]

    # --- exchange contract (statsys._thermal.ExchangeSubstrate) ---
    def neighbor_indices(self, nd: int) -> np.ndarray:
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        return self._nbr_idx[lo:hi]

    def _row_alignment(self, a: np.ndarray, b: np.ndarray) -> float:
        """``C(a, b) = Σ_{s1,s2} M[s1,s2] δ(a^{s1}, b^{s2})`` — the pair
        alignment of two site rows (symmetric in its arguments for
        symmetric M)."""
        eq = a[:, None] == b[None, :]
        return float((self.interaction_matrix * eq).sum())

    def _site_delta_E_row(self, nd: int, new_row: np.ndarray) -> float:
        """ΔE of replacing site *nd*'s ENTIRE row (all components at once):
        ``Σ_j J_nd,j [C(σ_nd, σ_j) − C(new_row, σ_j)]``."""
        lo, hi = self._nbr_ptr[nd], self._nbr_ptr[nd + 1]
        J = self._nbr_J[lo:hi]
        nbr_s = self.s[self._nbr_idx[lo:hi]]
        M = self.interaction_matrix
        cur = self.s[nd]
        # alignment of each neighbor row with cur / new_row
        eq_cur = cur[None, :, None] == nbr_s[:, None, :]
        eq_new = new_row[None, :, None] == nbr_s[:, None, :]
        c_cur = (eq_cur * M[None, :, :]).sum(axis=(1, 2))
        c_new = (eq_new * M[None, :, :]).sum(axis=(1, 2))
        return float((J * (c_cur - c_new)).sum())

    def swap_delta_E(self, u: int, v: int) -> float:
        """ΔE of swapping the FULL site rows at *u*, *v*: the two row ΔE's
        each miscount the shared bond (each replaces one endpoint by the
        other's row, producing self-alignment terms), but a swap leaves the
        pair alignment ``C(σ_u, σ_v)`` invariant — hence the
        ``+J_uv [C(σ_u,σ_u) + C(σ_v,σ_v) − 2 C(σ_u,σ_v)]`` correction (for
        k = 1, M = 1 this is the Potts ``+2 J_uv``)."""
        dE = self._site_delta_E_row(u, self.s[v]) + self._site_delta_E_row(
            v, self.s[u]
        )
        lo, hi = self._nbr_ptr[u], self._nbr_ptr[u + 1]
        row = self._nbr_idx[lo:hi]
        pos = np.flatnonzero(row == v)
        J_uv = float(self._nbr_J[lo:hi][pos[0]]) if len(pos) else 0.0
        corr = (
            self._row_alignment(self.s[u], self.s[u])
            + self._row_alignment(self.s[v], self.s[v])
            - 2.0 * self._row_alignment(self.s[u], self.s[v])
        )
        return dE + J_uv * corr

    def commit_swap(self, u: int, v: int) -> None:
        su = self.s[u].copy()
        self.s[u] = self.s[v]
        self.s[v] = su

    # ------------------------------------------------------------------
    # Energy (same J and M as delta_E, by construction)
    # ------------------------------------------------------------------
    def total_energy(self, s: np.ndarray | None = None) -> float:
        """Extensive ``H = −Σ_{(i,j)} Σ_{s1,s2} J_ij M[s1,s2]
        δ(σ_i^{s1}, σ_j^{s2})`` (each undirected edge once)."""
        self._ensure_couplings()
        state = self.s if s is None else np.asarray(s)
        eq = state[self._edge_u][:, :, None] == state[self._edge_v][:, None, :]
        pair_align = (eq * self.interaction_matrix[None, :, :]).sum(axis=(1, 2))
        return float(-(self._edge_J * pair_align).sum())

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
        self.observables[MULTISPEC_OBS_ENERGY].append(
            self._e_running / float(self.N)
        )

    def _record_magn(self) -> None:
        self.observables[MULTISPEC_OBS_MAGN].append(self.order_parameter())

    def _record_snapshot(self) -> None:
        snap = self.observables[MULTISPEC_OBS_SNAPSHOTS]
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
        row_len = self.species * self.N
        nbytes = n_rec_est * row_len * np.dtype(np.int32).itemsize
        if nbytes > MULTISPEC_SNAPSHOT_MAX_BYTES:
            raise ValueError(
                f"snapshot trajectory would be ~{nbytes / 2**20:.0f} MiB "
                f"({n_rec_est} snapshots x kN={row_len}, int32), exceeding "
                f"MULTISPEC_SNAPSHOT_MAX_BYTES="
                f"{MULTISPEC_SNAPSHOT_MAX_BYTES / 2**20:.0f} MiB. Increase "
                f"snapshot_every (currently {self.snapshot_every}), reduce "
                f"steps, or raise the cap in MultiSpeciesModel/defaults.py."
            )

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C: one directory per run)
    # ------------------------------------------------------------------
    def _physics_tokens(self) -> list:
        qs = list(self._q_per_species)
        q_val = qs[0] if len(set(qs)) == 1 else tuple(qs)
        M = self.interaction_matrix
        identity = bool(np.allclose(M, np.eye(self.species)))
        return [
            ("k", int(self.species)),
            ("q", q_val),
            ("T", float(self.T)),
            ("M", None if identity else "custom"),
        ]

    def _axis_tokens(self) -> list:
        """Dynamics axes (rule/move/upd/ord/tf) — leaf hook."""
        return []

    def _name_tokens(self) -> list:
        snap_on = MULTISPEC_OBS_SNAPSHOTS in self._selected_obs
        return [
            ("p", float(self.sg.pflip)),
            *self._physics_tokens(),
            ("ns", int(self.steps)),
            (
                "se",
                int(self.snapshot_every) if snap_on else None,
                MULTISPEC_SNAPSHOT_EVERY_DEFAULT,
            ),
            *self._axis_tokens(),
            ("cn", self.coupling_norm, MULTISPEC_COUPLING_NORM_DEFAULT),
            ("ic", self.ic, "uniform"),
            ("lang", self._lang_token()),
            ("s", int(self.seed)),
        ]

    def _begin_outputs(self) -> None:
        snap = self.observables[MULTISPEC_OBS_SNAPSHOTS]
        snap.reset_stream()
        self._snap_seen = 0
        if not self.savedisk or not self._is_py_runlang():
            return
        rundir = self._run_output_dir()
        rundir.mkdir(parents=True, exist_ok=True)
        self._write_cfg_sidecar(rundir)
        if MULTISPEC_OBS_ENERGY in self._selected_obs:
            self.observables[MULTISPEC_OBS_ENERGY].set_path(
                rundir / f"{MULTISPEC_ENERGY_FBASE}{BIN}"
            )
        if MULTISPEC_OBS_MAGN in self._selected_obs:
            self.observables[MULTISPEC_OBS_MAGN].set_path(
                rundir / f"{MULTISPEC_MAGN_FBASE}{BIN}"
            )
        if MULTISPEC_OBS_SNAPSHOTS in self._selected_obs:
            snap.set_path(rundir / f"{MULTISPEC_SNAPSHOTS_FBASE}{BIN}")
            self._guard_snapshot_size()
            snap.open_stream()

    def _persist_observables(self) -> None:
        if not self.savedisk or not self._is_py_runlang():
            return
        snap = self.observables[MULTISPEC_OBS_SNAPSHOTS]
        snap.close_stream()
        if (
            MULTISPEC_OBS_SNAPSHOTS in self._selected_obs
            and snap.path is not None
        ):
            snap.finalize_batch(subsample=self.snapshot_every)
        if MULTISPEC_OBS_ENERGY in self._selected_obs:
            self.observables[MULTISPEC_OBS_ENERGY].persist()
        if MULTISPEC_OBS_MAGN in self._selected_obs:
            self.observables[MULTISPEC_OBS_MAGN].persist()

    # ------------------------------------------------------------------
    # Engine plumbing (schemes provide the engine; the base runs it)
    # ------------------------------------------------------------------
    def _make_engine(self):
        raise NotImplementedError(
            "Scheme classes (MultiSpeciesMetropolis, ...) must implement "
            "_make_engine()."
        )

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """Python-backend sampling loop: record at sweep times 0..steps-1,
        advance one compiled sweep, track the energy incrementally from the
        sweep's ΔE (no per-sweep full recompute)."""
        engine = self._make_engine()
        self._engine = engine
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
            "The C-subprocess backend stays on the legacy MultiSpeciesModel "
            "class (frozen interface, decision D-B5); use runlang='py'/'pb' "
            "or the legacy MultiSpeciesModel."
        )
