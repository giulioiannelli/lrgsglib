"""RRG_TransCluster: topological percolation on random regular graphs.

Glue between the command line and the graph-agnostic kernels:

- ``resolve_run`` validates the parsed arguments and resolves them into an
  ``RRGRun``: disorder family, grid, tile layout, sampler, solver, run
  directory and identity. Invalid input raises ``ValueError``.
- ``run_transcluster_rrg`` publishes ``cfg.json`` once, then runs one tile
  (``-tid t``, a SLURM array task) or all of them (``-tid -1``).

Each rep draws its own topology (``kernels/RRG.draw_topology``) and sweeps
its tile's control values with ``kernels/TransClusterCore.sweep_rep``.

Run directory::

    $LRGSG_DATA/<wd>/phtra/N=<N>_k=<k>/<tokens>/{cfg.json, tiles/}

The tokens hold the identity (``mu``, ``tol``, ``cell``, ``ew``, ``ens``,
``smp``, ``ge`` for native samplers, ``tag``, ``lang``, ``s``) with
defaults elided. Grid, reps and tiling are not part of the identity, so
more grids and reps can be added to a run later.

Exit codes: ``TC_EXIT_OK``, ``TC_EXIT_CFG_MISMATCH`` (the run directory
holds another run) and ``TC_EXIT_PREEMPTED`` (signalled; resubmit).
Argument errors exit through ``argparse`` (code 2). Unreadable tile files
are set aside and recomputed (``kernels/TransClusterTiles.py``).
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray

from lrgsglib.config.const import PATHDATA, PATHNPHTR
from lrgsglib.config.progargs.defs.generic import DEFAULT_OUTSFFX
from lrgsglib.config.progargs.defs.RandomRegular import (
    DEFAULT_RRG_ENGINE,
    RRG_SAMPLER_PAIRING,
)
from lrgsglib.config.progargs.defs.TransCluster import (
    DEFAULT_TC_CELL,
    DEFAULT_TC_EDGE_WEIGHT,
    DEFAULT_TC_EIG_TOL,
    DEFAULT_TC_ENSEMBLE,
    DEFAULT_TC_MU,
    TC_ALL_TILES,
    TC_CELL_FAMILY,
    TC_CTRL_RANGE,
    TC_EDGE_WEIGHT_NORMAL,
    TC_EIG_TOL_MAX,
    TC_EXIT_CFG_MISMATCH,
    TC_EXIT_OK,
    TC_EXIT_PREEMPTED,
    TC_FAMILY_GAUSS,
    TC_LANG,
    TC_RUNDIR_NAME_MAX,
    TC_SEED_SCHEMA,
    TC_SEED_SPAN,
    TC_TAG_FORBIDDEN,
    TC_TILE_SCHEMA,
)
from lrgsglib.graphs import kRegularGraph
from lrgsglib.statsys._naming import (
    cfg_sidecar_bytes,
    publish_cfg_sidecar,
    run_dirname,
)

from .RRG import (
    draw_topology,
    identity_engine,
    require_engine,
    resolve_sampler,
    rrg_graph_key,
    validate_rrg,
)
from .TransClusterCore import Family, SolverConfig, check_ctrl, sweep_rep
from .TransClusterTiles import (
    Layout,
    Sweep,
    TileHeader,
    TileMismatch,
    decode_tile,
    git_state,
    grid_hash,
    make_layout,
    resolve_grid,
    run_tile,
    stop_signals,
    tile_paths,
    tile_stem,
)

__all__ = [
    "RRGRun",
    "family_from_args",
    "grid_from_args",
    "identity_blocks",
    "resolve_run",
    "rrg_rundir",
    "rrg_sweep",
    "run_transcluster_rrg",
    "tile_header",
    "tile_ids",
]


@dataclass(frozen=True)
class RRGRun:
    """Everything a task needs, resolved and validated."""

    n_nodes: int
    k: int
    family: Family
    grid: NDArray[np.float64]
    layout: Layout
    tile_ids: tuple[int, ...]
    seed: int
    engine: str
    method: str
    solver: SolverConfig
    gap_every: int
    rundir: Path
    identity: dict[str, dict[str, Any]]  # the graph / model / run blocks
    argv: tuple[str, ...]
    verbose: bool

    @property
    def identity_hash(self) -> str:
        """sha256 of the run's ``cfg.json``; every tile carries it."""
        return hashlib.sha256(cfg_sidecar_bytes(**self.identity)).hexdigest()


# ---------------------------------------------------------------------------
# Arguments -> run
# ---------------------------------------------------------------------------
def family_from_args(args: argparse.Namespace) -> Family:
    """Disorder family of ``-c`` / ``-ew``; knobs of other families raise."""
    if args.edge_weight == TC_EDGE_WEIGHT_NORMAL:
        if args.cell_type != DEFAULT_TC_CELL:
            raise ValueError("-c applies to -ew flip only")
        if args.ensemble != DEFAULT_TC_ENSEMBLE:
            raise ValueError("--ensemble applies to -ew flip only")
        return Family(TC_FAMILY_GAUSS, mu=args.mu)
    if args.mu != DEFAULT_TC_MU:
        raise ValueError("--mu applies to -ew normal only")
    return Family(TC_CELL_FAMILY[args.cell_type], ensemble=args.ensemble)


def grid_from_args(
    args: argparse.Namespace, family: Family
) -> NDArray[np.float64]:
    """Control grid of the one grid flag given; it must fit the family."""
    gauss = family.kind == TC_FAMILY_GAUSS
    p_values = args.p_list or args.p_linsp  # nargs="+" lists are never empty
    sigmas = args.sigma_list or args.sigma_linsp
    if gauss and sigmas is None:
        raise ValueError("-ew normal takes --sigma-list or --sigma-linsp")
    if not gauss and p_values is None:
        raise ValueError("-ew flip takes --p-list or --p-linsp")
    low, high = TC_CTRL_RANGE[family.kind]
    grid = resolve_grid(sigmas if gauss else p_values, low, high)
    for ctrl in grid:
        check_ctrl(family, ctrl)  # sigma = 0 at mu = 0 has no couplings
    return grid


def tile_ids(layout: Layout, tile_id: int) -> tuple[int, ...]:
    """The tiles of this call: one, or all of them for ``TC_ALL_TILES``."""
    if tile_id == TC_ALL_TILES:
        return tuple(range(layout.n_tiles))
    decode_tile(layout, tile_id)  # raises if out of range
    return (tile_id,)


def identity_blocks(
    n_nodes: int,
    k: int,
    family: Family,
    method: str,
    engine: str,
    seed: int,
    solver: SolverConfig,
) -> dict[str, dict[str, Any]]:
    """The ``cfg.json`` identity: what the ensemble and the numbers depend on.

    Grid, reps, tiling and ``--gap-every`` stay out, so a run can grow.
    """
    gauss = family.kind == TC_FAMILY_GAUSS
    graph = {
        "type": kRegularGraph.__name__,
        "N": n_nodes,
        "k": k,
        "method": method,
        "engine": identity_engine(method, engine),
        "seed": seed,
    }
    model = {
        "family": family.kind,
        "ensemble": None if gauss else family.ensemble,
        "mu": family.mu if gauss else None,
    }
    run = {
        "tile_schema": TC_TILE_SCHEMA,
        "seed_schema": TC_SEED_SCHEMA,
        "lang": TC_LANG,
        "solver": dataclasses.asdict(solver),
    }
    return {"graph": graph, "model": model, "run": run}


def rrg_rundir(args: argparse.Namespace, family: Family, method: str) -> Path:
    """Run directory of this identity (see the module docstring). A name of
    more than ``TC_RUNDIR_NAME_MAX`` bytes raises ``ValueError``."""
    gauss = family.kind == TC_FAMILY_GAUSS
    tokens = [
        ("mu", family.mu if gauss else None, DEFAULT_TC_MU),
        ("tol", args.eig_tol, DEFAULT_TC_EIG_TOL),
        ("cell", None if gauss else args.cell_type, DEFAULT_TC_CELL),
        ("ew", args.edge_weight, DEFAULT_TC_EDGE_WEIGHT),
        ("ens", None if gauss else family.ensemble, DEFAULT_TC_ENSEMBLE),
        ("smp", method, RRG_SAMPLER_PAIRING),
        ("ge", identity_engine(method, args.graph_engine), DEFAULT_RRG_ENGINE),
        ("tag", args.out_suffix, DEFAULT_OUTSFFX),
        ("lang", TC_LANG),
        ("s", args.seed),
    ]
    name = run_dirname(tokens)
    size = len(name.encode())
    if size > TC_RUNDIR_NAME_MAX:
        raise ValueError(
            f"the run directory name takes {size} bytes, more than "
            f"{TC_RUNDIR_NAME_MAX}: shorten -os"
        )
    shape = f"N={args.N}_k={args.k}"
    return PATHDATA / args.workdir / PATHNPHTR / shape / name


def _check_scalars(args: argparse.Namespace) -> None:
    """Range checks the parser cannot express."""
    if args.seed < 0:
        raise ValueError("--seed must be >= 0")
    if not 0 <= args.eig_tol <= TC_EIG_TOL_MAX:  # NaN fails as well
        raise ValueError(f"--eig-tol must lie in [0, {TC_EIG_TOL_MAX:g}]")
    if args.gap_every < 0:
        raise ValueError("--gap-every must be >= 0")
    if args.rep_start + args.number_of_averages > TC_SEED_SPAN:
        raise ValueError(f"rep indices must stay below {TC_SEED_SPAN}")
    if any(char in args.out_suffix for char in TC_TAG_FORBIDDEN):
        raise ValueError(f"-os cannot contain any of {TC_TAG_FORBIDDEN!r}")


def resolve_run(args: argparse.Namespace, argv: Sequence[str]) -> RRGRun:
    """Validate ``args`` and resolve them into an ``RRGRun``."""
    validate_rrg(args.N, args.k)
    _check_scalars(args)
    # -0.0 -> 0.0 in the float identity values: one run dir, one cfg.json
    zeros = {"eig_tol": args.eig_tol + 0.0, "mu": args.mu + 0.0}
    args = argparse.Namespace(**{**vars(args), **zeros})
    family = family_from_args(args)
    grid = grid_from_args(args, family)
    layout = make_layout(
        grid.size,
        args.rep_start,
        args.number_of_averages,
        args.ctrl_chunk,
        args.rep_chunk,
    )
    require_engine(args.graph_engine)
    method = resolve_sampler(args.sampler, args.k, args.graph_engine)
    solver = SolverConfig(tol=args.eig_tol)
    identity = identity_blocks(
        args.N, args.k, family, method, args.graph_engine, args.seed, solver
    )
    cfg_sidecar_bytes(**identity)  # NaN or inf: exit 2 here, not 1 at publish
    return RRGRun(
        n_nodes=args.N,
        k=args.k,
        family=family,
        grid=grid,
        layout=layout,
        tile_ids=tile_ids(layout, args.tile_id),
        seed=args.seed,
        engine=args.graph_engine,
        method=method,
        solver=solver,
        gap_every=args.gap_every,
        rundir=rrg_rundir(args, family, method),
        identity=identity,
        argv=tuple(argv),
        verbose=args.verbose,
    )


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def rrg_sweep(run: RRGRun) -> Sweep:
    """``sweep(rep, ictrls)``: draw rep's topology, then measure the units."""
    key = rrg_graph_key(run.n_nodes, run.k)

    def sweep(rep: int, ictrls: Sequence[int]):
        topo = draw_topology(
            run.n_nodes,
            run.k,
            root=run.seed,
            rep=rep,
            engine=run.engine,
            method=run.method,
        )
        base = {"attempt": topo.attempt, "topo_fp": topo.fingerprint}
        return sweep_rep(
            topo.edges,
            run.n_nodes,
            run.family,
            run.grid,
            ictrls,
            root=run.seed,
            graph_key=key,
            rep=rep,
            gap_every=run.gap_every,
            solver=run.solver,
            base_fields=base,
        )

    return sweep


def tile_header(run: RRGRun) -> TileHeader:
    """What every tile of ``run`` carries besides its units."""
    return TileHeader(run.identity_hash, run.grid, run.gap_every, run.argv)


def run_transcluster_rrg(run: RRGRun) -> int:
    """Publish ``cfg.json`` once, then run the tiles; return an exit code.

    One stop flag covers every tile: a signal during a tile's last unit
    still ends the run before the next tile.
    """
    git_state()  # the provenance describes the code imported at start
    with stop_signals() as stop:
        try:
            publish_cfg_sidecar(run.rundir, **run.identity)
        except FileExistsError as err:
            print(f"error: {err}", file=sys.stderr)
            return TC_EXIT_CFG_MISMATCH
        header, ghash = tile_header(run), grid_hash(run.grid)
        sweep = rrg_sweep(run)
        if run.verbose:
            print(f"run dir {run.rundir} ({len(run.tile_ids)} tiles)")
        for tile_id in run.tile_ids:
            if stop.signum is not None:
                return TC_EXIT_PREEMPTED
            tile, start = decode_tile(run.layout, tile_id), time.perf_counter()
            paths = tile_paths(run.rundir, ghash, tile)
            try:
                code = run_tile(sweep, tile, paths, header, stop)
            except TileMismatch as err:
                print(f"error: {err}", file=sys.stderr)
                return TC_EXIT_CFG_MISMATCH
            if code != TC_EXIT_OK:
                return code
            if run.verbose:
                elapsed = time.perf_counter() - start
                print(
                    f"tile {tile_id} {tile_stem(ghash, tile)}: {elapsed:.1f} s"
                )
    return TC_EXIT_OK
