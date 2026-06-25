"""Signed random-walker dynamics — single-walker simulator family.

Each walker is a particle moving on the nodes of a :class:`SignedGraph`.
Three rule variants are shipped (:class:`AbsorbingWalker`,
:class:`KillingWalker`, :class:`StickyWalker`); X-node behaviour is
controlled by a ``x_node_behavior`` kwarg (``'reflect'`` by default —
physically: unfrustrated X-nodes act as reflective barriers, not as
absorbers).

This module only provides the single-walker base class. Higher-level
experiments (two-walker overlap, disorder averaging) live in
:mod:`_overlap.py` and in the corresponding CLI program
``src/L2D_SignedRW.py``; the walker class is deliberately agnostic to
them, mirroring the separation between :class:`IsingDynamics` and
the downstream specific-heat / susceptibility analyses.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Tuple

import numpy as np

from ...config.const import (
    DEFAULT_SRW_COVERAGE_FRAC,
    DEFAULT_SRW_MAX_N_CROSS,
    DEFAULT_SRW_N_WALKERS,
    DEFAULT_SRW_SEED,
    DEFAULT_SRW_START,
    DEFAULT_SRW_X_NODE,
    SRW_START_PROTOCOLS,
    SRW_X_NODE_BEHAVIORS,
)
from ..DynSys import DynSys
from .._solver import SolverBackend
from .._solver_engine import get_solver
from ._kernel import _kill_masks, run_walker, signed_lattice_tables
from .defaults import SRW_WALKER_SOLVER_NAME

if TYPE_CHECKING:
    from ...graphs.protocols import SignedGraphProtocol as SignedGraph


def _load_srw_native():
    """Lazy-import the compiled _srw_native extension.

    Returns None if the shared library is not built.
    """
    try:
        from .ccore import _srw_native  # type: ignore
        return _srw_native
    except ImportError:
        try:
            # When the .so lives next to ccore/__init__.py but is not
            # exported there, add the directory to sys.path manually.
            import importlib.util
            from pathlib import Path
            so_dir = Path(__file__).resolve().parent / "ccore"
            candidates = list(so_dir.glob("_srw_native*.so"))
            if not candidates:
                return None
            spec = importlib.util.spec_from_file_location(
                "_srw_native", candidates[0]
            )
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
        except Exception:
            return None


class SignedWalker(DynSys):
    """Abstract base for the signed-walker family.

    Subclasses set :attr:`rule` to one of ``{'absorb', 'kill', 'sticky'}``
    and otherwise inherit everything.

    Parameters
    ----------
    sg : SignedGraph
        Any object satisfying :class:`SignedGraphProtocol` (NX or GT).
    n_walkers : int
    seed : int
    start : {'random', 'fixed', 'center'}
        Starting-position protocol. ``'fixed'`` requires
        ``start_kwargs={'node': int}`` or ``{'nodes': array-like}``.
        ``'center'`` is lattice-aware: uses index ``(s//2)*(s+1)`` where
        ``s = round(N**0.5)``.
    start_kwargs : dict, optional
    coverage_stop : float
        Terminate a walker once ``unique_visits / N >= coverage_stop``.
    x_node_behavior : {'reflect', 'absorb'}
        Fate of unfrustrated negative edges (edges incident to X-nodes).
        ``'reflect'`` (library default) treats them as reflective
        barriers; ``'absorb'`` matches the notebook's simpler rule
        (any negative edge crossing kills).
    max_n_cross : int
        Sticky-rule only — death threshold for the per-walker frustrated
        crossings counter. Ignored by absorb/kill rules.
    store_trajectory : bool
        If True, ``.trajectory`` holds the full ``(T+1, n_walkers)``
        history. Memory: ``8 * n_walkers * T`` bytes.
    store_per_walker_visits : bool
        If True, ``.visits_per_walker`` holds an ``(n_walkers, N)``
        int32 array. Memory: ``4 * n_walkers * N`` bytes.
    runlang : str
        Only ``'py'`` is implemented in Phase 1; ``'pb_*'`` / ``'cu_*'``
        / ``'c_*'`` are reserved for Phases 2-3.
    """

    rule: str = ""  # overridden by subclasses
    dyn_UVclass: str = "signed_walker"

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        n_walkers: int = DEFAULT_SRW_N_WALKERS,
        seed: int = DEFAULT_SRW_SEED,
        start: str = DEFAULT_SRW_START,
        start_kwargs: Optional[dict] = None,
        coverage_stop: float = DEFAULT_SRW_COVERAGE_FRAC,
        x_node_behavior: str = DEFAULT_SRW_X_NODE,
        max_n_cross: int = DEFAULT_SRW_MAX_N_CROSS,
        store_trajectory: bool = False,
        store_per_walker_visits: bool = False,
        runlang: str = "py",
        **kwargs: Any,
    ) -> None:
        if self.rule == "":
            raise TypeError(
                f"{type(self).__name__} is abstract; instantiate a concrete "
                f"subclass (AbsorbingWalker, KillingWalker, StickyWalker)."
            )
        if start not in SRW_START_PROTOCOLS:
            raise ValueError(
                f"start must be one of {SRW_START_PROTOCOLS}, got {start!r}."
            )
        if x_node_behavior not in SRW_X_NODE_BEHAVIORS:
            raise ValueError(
                f"x_node_behavior must be one of {SRW_X_NODE_BEHAVIORS}, "
                f"got {x_node_behavior!r}."
            )
        if not 0.0 < coverage_stop <= 1.0:
            raise ValueError(
                f"coverage_stop must be in (0, 1], got {coverage_stop!r}."
            )
        if n_walkers < 1:
            raise ValueError(f"n_walkers must be >= 1, got {n_walkers!r}.")

        self.n_walkers = int(n_walkers)
        self.start = start
        self.start_kwargs = dict(start_kwargs) if start_kwargs else {}
        self.coverage_stop = float(coverage_stop)
        self.x_node_behavior = x_node_behavior
        self.max_n_cross = int(max_n_cross)
        self.store_trajectory = bool(store_trajectory)
        self.store_per_walker_visits = bool(store_per_walker_visits)

        dynpath = getattr(sg, "path_srw", None)
        # walker time is not N-sweep time — pass steps=1 as placeholder so the
        # DynSys time-control helper is happy; real per-walker stop steps are
        # recorded in ``self.stop_step`` after ``run()``.
        super().__init__(sg, runlang=runlang, steps=1, seed=seed,
                         dynpath=dynpath, **kwargs)

    # ------------------------------------------------------------------
    # DynSys abstract API (minimal stubs; walker uses vectorised kernel)
    # ------------------------------------------------------------------
    _state_dtype = np.int64  # per-walker node index

    @property
    def _state_shape(self) -> Tuple[int, ...]:
        return (self.n_walkers,)

    def init_state(self, custom: Any = None) -> None:
        """Populate ``self.s`` with initial walker positions."""
        self.s = self._resolve_start_positions().copy()

    def export_state(self) -> None:
        raise NotImplementedError(
            "SignedWalker has no C-backend state export yet (Phase 1)."
        )

    def ds1step(self, nd: int) -> None:
        raise NotImplementedError(
            "SignedWalker uses a vectorised kernel; ds1step is not applicable."
        )

    def run_c(self, **kw: Any) -> None:
        raise NotImplementedError(
            "SignedWalker C backend is deferred to Phase 2/3."
        )

    # ------------------------------------------------------------------
    # Start-position protocols
    # ------------------------------------------------------------------
    def _resolve_start_positions(self) -> np.ndarray:
        """Compute the starting-node index array according to ``self.start``."""
        N = int(self.sg.N)
        rng = np.random.default_rng(self.seed)
        if self.start == 'random':
            return rng.integers(0, N, size=self.n_walkers, dtype=np.int64)
        if self.start == 'fixed':
            if 'nodes' in self.start_kwargs:
                nodes = np.asarray(self.start_kwargs['nodes'], dtype=np.int64)
                if nodes.shape == (self.n_walkers,):
                    return nodes
                if nodes.ndim == 0:
                    return np.full(self.n_walkers, int(nodes), dtype=np.int64)
                raise ValueError(
                    f"start_kwargs['nodes'] shape {nodes.shape} "
                    f"!= ({self.n_walkers},) or scalar."
                )
            if 'node' in self.start_kwargs:
                return np.full(self.n_walkers,
                               int(self.start_kwargs['node']),
                               dtype=np.int64)
            raise ValueError(
                "start='fixed' requires start_kwargs={'node': int} "
                "or start_kwargs={'nodes': array-like}."
            )
        if self.start == 'center':
            side = int(round(N ** 0.5))
            if side * side != N:
                raise ValueError(
                    f"start='center' requires a square lattice "
                    f"(N={N} is not a perfect square)."
                )
            center_idx = (side // 2) * (side + 1)
            return np.full(self.n_walkers, center_idx, dtype=np.int64)
        raise ValueError(f"unknown start protocol {self.start!r}.")

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------
    def _resolve_backend(self) -> SolverBackend:
        """Map the runlang code to a solver family (``py`` / ``pb_absorb``)."""
        if self.runlang.startswith('pb'):
            return SolverBackend.PB
        if self.runlang.startswith('py'):
            return SolverBackend.PY
        raise NotImplementedError(
            f"runlang={self.runlang!r} is not available yet; "
            f"use 'py' or 'pb_absorb'."
        )

    def run(self, **kw: Any) -> None:
        """Dispatch to the chosen backend through the shared solver registry.

        Supported runlang values:

        * ``'py'`` — reference numpy kernel (all three rules).
        * ``'pb_absorb'`` — pybind11 / C kernel for the absorb rule (requires
          ``rule='absorb'``). Statistical parity with the Python backend
          (different RNG).
        """
        backend = self._resolve_backend()
        solver = get_solver(SRW_WALKER_SOLVER_NAME, backend)
        solver.supports(self)
        solver.execute(self, **kw)

    def run_py(self, **kw: Any) -> None:
        """Run the Python kernel and populate per-walker observables."""
        self.init_state()
        result = run_walker(
            self.sg, self.rule,
            n_walkers=self.n_walkers,
            start_positions=self.s,
            seed=self.seed,
            coverage_stop=self.coverage_stop,
            x_node_behavior=self.x_node_behavior,
            max_n_cross=self.max_n_cross,
            store_trajectory=self.store_trajectory,
            store_per_walker_visits=self.store_per_walker_visits,
        )
        self._bind_result(result)

    def _run_pybind(self, **kw: Any) -> None:
        """Run the C absorb-kernel via the pybind11 extension."""
        if self.runlang != 'pb_absorb':
            raise NotImplementedError(
                f"Phase 2 ships only 'pb_absorb' (not {self.runlang!r}); "
                f"kill/sticky C kernels are a follow-up."
            )
        if self.rule != 'absorb':
            raise RuntimeError(
                f"runlang='pb_absorb' requires rule='absorb'; "
                f"got rule={self.rule!r}."
            )
        native = _load_srw_native()
        if native is None:
            raise RuntimeError(
                "_srw_native extension is not built; run `make` under "
                "lrgsglib/src/lrgsglib/statsys/SignedRW/ccore/native/."
            )
        self.init_state()
        # Build CSR-style mask arrays from the lattice's (N, z) tables.
        nbrs, neg, frust = signed_lattice_tables(self.sg)
        kill_mask, reflect_mask = _kill_masks(neg, frust, self.x_node_behavior)
        N, z = nbrs.shape
        stop_thresh = max(1, int(self.coverage_stop * N))
        neigh_indices = np.ascontiguousarray(nbrs.ravel(), dtype=np.int64)
        neigh_ptr = np.arange(0, (N + 1) * z, z, dtype=np.int64)
        kill_flat = np.ascontiguousarray(kill_mask.ravel(), dtype=np.uint8)
        refl_flat = np.ascontiguousarray(reflect_mask.ravel(), dtype=np.uint8)
        start_positions = np.ascontiguousarray(self.s, dtype=np.int64)
        # walkers_per_trial=1 recovers a per-walker bitmap via the trial
        # path; 0 disables the optional output entirely.
        wpt = 1 if self.store_per_walker_visits else 0
        (
            visits_agg, unique_visits, stop_step, stop_reason,
            final_position, trial_bubbles,
        ) = native.absorb_walker_sampling(
            N,
            neigh_indices,
            neigh_ptr,
            kill_flat,
            refl_flat,
            start_positions,
            int(stop_thresh),
            int(self.seed),
            wpt,
        )
        Ne = int(self.sg.Ne) if hasattr(self.sg, 'Ne') else int(self.sg.num_edges)
        # trial_bubbles is (n_walkers, N) when walkers_per_trial=1, so it
        # coincides with the per-walker visited bitmap.
        vpw = trial_bubbles if self.store_per_walker_visits else None
        result = dict(
            unique_visits=unique_visits.astype(np.int64, copy=False),
            unique_frac=unique_visits / N,
            stop_step=stop_step.astype(np.int64, copy=False),
            stop_reason=stop_reason.astype(np.int8, copy=False),
            final_position=final_position.astype(np.int64, copy=False),
            visits_agg=visits_agg.astype(np.int64, copy=False),
            visits_per_walker=vpw,
            history=None,
            neg_frac=float(neg.sum()) / (2 * Ne),
            frust_frac=float(frust.sum()) / (2 * Ne),
        )
        self._bind_result(result)

    def _bind_result(self, result: dict) -> None:
        """Copy a kernel-result dict onto the walker instance."""
        self.unique_visits = result['unique_visits']
        self.unique_frac = result['unique_frac']
        self.stop_step = result['stop_step']
        self.stop_reason = result['stop_reason']
        self.final_position = result['final_position']
        self.visits_agg = result['visits_agg']
        self.visits_per_walker = result['visits_per_walker']
        self.trajectory = result['history']
        self.neg_frac = result['neg_frac']
        self.frust_frac = result['frust_frac']
        self.result = result

    # ------------------------------------------------------------------
    # Convenience summary
    # ------------------------------------------------------------------
    def summary(self) -> dict:
        """Return a picklable dict of the per-walker observables."""
        if not hasattr(self, 'unique_visits'):
            raise RuntimeError("run() must be called before summary().")
        return dict(
            rule=self.rule,
            n_walkers=self.n_walkers,
            coverage_stop=self.coverage_stop,
            x_node_behavior=self.x_node_behavior,
            start=self.start,
            unique_visits=self.unique_visits,
            unique_frac=self.unique_frac,
            stop_step=self.stop_step,
            stop_reason=self.stop_reason,
            final_position=self.final_position,
            visits_agg=self.visits_agg,
            neg_frac=self.neg_frac,
            frust_frac=self.frust_frac,
        )


class AbsorbingWalker(SignedWalker):
    """D3 — walker dies crossing a killing edge; reflects at unfrustrated negatives."""
    rule = 'absorb'
    dyn_UVclass = "srw_absorb"


class KillingWalker(SignedWalker):
    """D2 — walker dies with prob ``1/n_kill`` at any node with incident killing edges."""
    rule = 'kill'
    dyn_UVclass = "srw_kill"


class StickyWalker(SignedWalker):
    """D1 — walker's move prob halves per frustrated crossing; dies at ``max_n_cross``."""
    rule = 'sticky'
    dyn_UVclass = "srw_sticky"


__all__ = [
    'SignedWalker',
    'AbsorbingWalker',
    'KillingWalker',
    'StickyWalker',
]
