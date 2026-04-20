"""Kernel for ``L2D_SignedRW.py``.

Builds a 2D signed lattice (engine-aware via the standard ``-ge`` flag),
runs one (single mode) or two (pair mode) walker ensembles, writes per-
walker observables to the lattice's ``path_srw`` directory, and — in
pair mode — additionally writes a JSON summary of every overlap metric.

The walker class does NOT perform any overlap computation; overlap is
strictly an analysis step here, mirroring how the Ising program reports
energies/magnetizations while susceptibilities / specific heat are
computed downstream.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Any

import numpy as np

from lrgsglib.statsys.SignedRW import (
    AbsorbingWalker,
    KillingWalker,
    StickyWalker,
)
from lrgsglib.statsys.SignedRW._overlap import overlap_all

from .L2D import prepare_lattice as _prepare_lattice


_RULE_CLASSES = {
    'absorb': AbsorbingWalker,
    'kill': KillingWalker,
    'sticky': StickyWalker,
}


def _build_walker(args: Any, lattice: Any, *, seed: int,
                  start: str, start_node: int | None) -> Any:
    """Instantiate the walker class selected by ``args.rule``."""
    Cls = _RULE_CLASSES[args.rule]
    start_kwargs = {'node': start_node} if start == 'fixed' else None
    return Cls(
        lattice,
        n_walkers=args.n_walkers,
        seed=seed,
        start=start,
        start_kwargs=start_kwargs,
        coverage_stop=args.coverage,
        x_node_behavior=args.x_node,
        max_n_cross=args.max_n_cross,
        store_trajectory=args.store_trajectory,
        store_per_walker_visits=args.store_per_walker_visits,
        runlang=args.runlang,
    )


def _output_dir(lattice: Any, rule: str) -> Path:
    """Return the rule-specific subdirectory under ``path_srw``."""
    base = Path(lattice.path_srw) / rule
    base.mkdir(parents=True, exist_ok=True)
    return base


def _run_stem(args: Any, *, tag: str, seed: int) -> str:
    """Build a deterministic filename stem for a run."""
    from lrgsglib.config.funcs import peq_fstr
    # peq_fstr already returns 'p=<value>'; use it directly.
    return f"{tag}_{peq_fstr(args.p)}_L={args.L}_seed={seed}"


def _save_walker_pickle(walker: Any, out_dir: Path, stem: str,
                        verbose: bool = False) -> Path:
    """Persist the walker's summary as a pickle under ``out_dir``."""
    payload = walker.summary()
    payload['trajectory'] = walker.trajectory  # may be None
    payload['visits_per_walker'] = walker.visits_per_walker  # may be None
    path = out_dir / f"{stem}.pkl"
    with path.open('wb') as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    if verbose:
        print(f"  wrote {path}")
    return path


def _save_overlap_json(c_a: np.ndarray, c_b: np.ndarray,
                       out_dir: Path, stem: str,
                       verbose: bool = False) -> Path:
    """Compute every overlap flavour and dump them as JSON."""
    payload = overlap_all(c_a, c_b)
    path = out_dir / f"{stem}_overlap.json"
    with path.open('w') as fh:
        json.dump(payload, fh, indent=2)
    if verbose:
        print(f"  wrote {path}  -- l2_raw={payload['l2_raw']:.4g}")
    return path


def run_simulation(args: Any) -> None:
    """Entry point wired from ``src/L2D_SignedRW.py``."""
    lattice = _prepare_lattice(args)
    out_dir = _output_dir(lattice, args.rule)
    verbose = bool(getattr(args, 'verbose', False))

    # Walker A — always run, in both modes
    walker_a = _build_walker(
        args, lattice,
        seed=args.seed_a, start=args.start_a, start_node=args.start_a_node,
    )
    walker_a.run()
    stem_a = _run_stem(args, tag='walkerA', seed=args.seed_a)
    _save_walker_pickle(walker_a, out_dir, stem_a, verbose=verbose)

    if args.mode == 'pair':
        walker_b = _build_walker(
            args, lattice,
            seed=args.seed_b, start=args.start_b, start_node=args.start_b_node,
        )
        walker_b.run()
        stem_b = _run_stem(args, tag='walkerB', seed=args.seed_b)
        _save_walker_pickle(walker_b, out_dir, stem_b, verbose=verbose)
        stem_pair = _run_stem(args, tag='pair', seed=args.seed_a)
        _save_overlap_json(
            walker_a.visits_agg, walker_b.visits_agg,
            out_dir, stem_pair, verbose=verbose,
        )

    if verbose:
        print(
            f"L2D_SignedRW done — rule={args.rule} mode={args.mode} "
            f"L={args.L} p={args.p:.4g} "
            f"x_node={args.x_node} n_walkers={args.n_walkers}"
        )
