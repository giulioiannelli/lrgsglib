"""Tiles, tile files and the resumable tile loop of TransCluster v2 programs.

A *unit* is one ``(rep, ictrl)`` pair (see ``TransClusterCore``). A *tile* is
the rectangle of units with reps in ``[r0, r1)`` and control indices in
``[c0, c1)``. Tile ids decode rep-major: id ``t`` covers rep block
``t // n_ctrl_tiles`` and control block ``t % n_ctrl_tiles``, so a throttled
array finishes whole low-statistics curves first. With one tile per task,
``ctrl_chunk = 1`` splits a run over the control value, ``rep_chunk = 1``
over the disorder, and both together give one unit per task.

Tile files live in ``<rundir>/tiles/``, with half-open ranges as in Python
slices::

    g=<grid hash>_c=<c0>-<c1>_r=<r0>-<r1>.npz        a complete tile
    g=<grid hash>_c=<c0>-<c1>_r=<r0>-<r1>.part.npz   checkpoint of a tile

A tile npz holds the record rows (``RECORD_DTYPE``, sorted by rep then
ictrl), the P(S) of each row in CSR form (``ps_ptr``, ``ps_size``,
``ps_count``), the grid, the tile bounds, the run identity (sha256 of
``cfg.json``), a completion flag and the provenance JSON. Records carry
their control value, so tiles of several grids can share a run directory.

``run_tile`` resumes from the partial and skips finished units. It writes a
checkpoint every ``TC_CKPT_SECONDS``, and once the caller's ``StopFlag``
holds SIGUSR1/SIGTERM it checkpoints and returns ``TC_EXIT_PREEMPTED``. A
tile file that cannot be read, or a partial of another run, is renamed to
``<name>.<id>.bad`` and the tile is recomputed: records only depend on the
seed. A complete tile of another run is an error. Complete tiles are
hard-linked into place (``publish_tile``), so a task never replaces the
complete tile that another task published while it ran.
"""

from __future__ import annotations

import functools
import hashlib
import io
import itertools
import json
import os
import platform
import secrets
import shlex
import signal
import socket
import subprocess
import sys
import time
import warnings
import zipfile
import zlib
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Optional, Sequence

import numpy as np
from numpy.lib.npyio import NpzFile
from numpy.typing import NDArray

from lrgsglib.config.const import EXACT_FLOAT_FMT
from lrgsglib.config.progargs.defs.TransCluster import (
    TC_BAD_EXT,
    TC_BAD_ID_BYTES,
    TC_CKPT_SECONDS,
    TC_EXIT_OK,
    TC_EXIT_PREEMPTED,
    TC_GIT_TIMEOUT,
    TC_GRID_HASH_BYTES,
    TC_PART_EXT,
    TC_PROVENANCE_MODULES,
    TC_SLURM_ENV_VARS,
    TC_STOP_SIGNALS,
    TC_TILE_EXT,
    TC_TILE_SCHEMA,
    TC_TILES_DIR,
)
from lrgsglib.statsys._io import atomic_write_bytes, create_once

from .TransClusterCore import PS, RECORD_DTYPE

__all__ = [
    "Layout",
    "StopFlag",
    "Sweep",
    "Tile",
    "TILE_KEYS",
    "TileBuffer",
    "TileError",
    "TileHeader",
    "TileMismatch",
    "TilePaths",
    "TileUnreadable",
    "all_tiles",
    "check_tile",
    "concat_int64",
    "decode_tile",
    "git_state",
    "grid_hash",
    "load_tile",
    "make_layout",
    "pending_units",
    "provenance",
    "publish_tile",
    "read_tile",
    "resolve_grid",
    "resume_buffer",
    "run_tile",
    "set_aside",
    "stop_signals",
    "tile_paths",
    "tile_stem",
    "write_npz",
    "write_tile",
]

#: ``sweep(rep, ictrls)`` yields the ``(record, P(S))`` of those units.
Sweep = Callable[[int, Sequence[int]], Iterable[tuple[np.ndarray, PS]]]


# ---------------------------------------------------------------------------
# Control grid
# ---------------------------------------------------------------------------
def resolve_grid(
    values: Iterable[float], low: float, high: float
) -> NDArray[np.float64]:
    """Sorted, deduplicated control grid, rounded to ``EXACT_FLOAT_FMT``.

    The rounding removes float noise (``linspace`` can give
    0.30000000000000004) and -0.0 becomes 0.0, so a point shared by two
    grids is the same float and its units give the same records.
    """
    grid = np.unique([float(format(v, EXACT_FLOAT_FMT)) + 0.0 for v in values])
    if grid.size == 0:
        raise ValueError("the control grid is empty (check the linspace num)")
    if not np.all(np.isfinite(grid)):
        raise ValueError("the control grid needs finite values")
    if grid[0] < low or grid[-1] > high:
        raise ValueError(f"control values must lie in [{low}, {high}]")
    return grid


def grid_hash(grid: NDArray[np.float64]) -> str:
    """Short hex digest of a resolved grid, used in tile names."""
    data = np.ascontiguousarray(grid, dtype="<f8").tobytes()
    return hashlib.blake2b(data, digest_size=TC_GRID_HASH_BYTES).hexdigest()


# ---------------------------------------------------------------------------
# Layout and tiles
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Layout:
    """How the units of a run split into tiles (chunk sizes resolved)."""

    n_ctrl: int
    rep_start: int
    n_reps: int
    ctrl_chunk: int
    rep_chunk: int

    @property
    def n_ctrl_tiles(self) -> int:
        return -(-self.n_ctrl // self.ctrl_chunk)

    @property
    def n_rep_tiles(self) -> int:
        return -(-self.n_reps // self.rep_chunk)

    @property
    def n_tiles(self) -> int:
        return self.n_ctrl_tiles * self.n_rep_tiles


def make_layout(
    n_ctrl: int, rep_start: int, n_reps: int, ctrl_chunk: int, rep_chunk: int
) -> Layout:
    """Layout of a run; a chunk of 0 (or larger than the axis) means all."""
    if n_ctrl < 1 or n_reps < 1 or rep_start < 0:
        raise ValueError("need >= 1 control value, >= 1 rep, rep start >= 0")
    if ctrl_chunk < 0 or rep_chunk < 0:
        raise ValueError("chunk sizes must be >= 0 (0 = all)")
    ctrl_chunk = min(ctrl_chunk or n_ctrl, n_ctrl)
    rep_chunk = min(rep_chunk or n_reps, n_reps)
    return Layout(n_ctrl, rep_start, n_reps, ctrl_chunk, rep_chunk)


@dataclass(frozen=True)
class Tile:
    """Units with rep in ``[r0, r1)`` and control index in ``[c0, c1)``."""

    c0: int
    c1: int
    r0: int
    r1: int

    @property
    def size(self) -> int:
        return (self.c1 - self.c0) * (self.r1 - self.r0)

    def units(self) -> set[tuple[int, int]]:
        """Every ``(rep, ictrl)`` of the tile."""
        reps, ctrls = range(self.r0, self.r1), range(self.c0, self.c1)
        return set(itertools.product(reps, ctrls))

    def bounds(self) -> NDArray[np.int64]:
        return np.array([self.c0, self.c1, self.r0, self.r1], dtype=np.int64)


def decode_tile(layout: Layout, tile_id: int) -> Tile:
    """Tile ``tile_id`` of ``layout``, decoded rep-major."""
    if not 0 <= tile_id < layout.n_tiles:
        raise ValueError(
            f"tile id must lie in [0, {layout.n_tiles}), got {tile_id}"
        )
    ir, ic = divmod(tile_id, layout.n_ctrl_tiles)
    c0 = ic * layout.ctrl_chunk
    r0 = layout.rep_start + ir * layout.rep_chunk
    c1 = min(c0 + layout.ctrl_chunk, layout.n_ctrl)
    r1 = min(r0 + layout.rep_chunk, layout.rep_start + layout.n_reps)
    return Tile(c0, c1, r0, r1)


def all_tiles(layout: Layout) -> list[Tile]:
    """Every tile of ``layout`` in id order; they cover each unit once."""
    tiles = [decode_tile(layout, t) for t in range(layout.n_tiles)]
    covered = sum(tile.size for tile in tiles)
    if covered != layout.n_ctrl * layout.n_reps:
        raise RuntimeError(f"tiles cover {covered} units, not all of them")
    return tiles


# ---------------------------------------------------------------------------
# Tile files
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TilePaths:
    """The complete tile file and its checkpoint."""

    final: Path
    part: Path


@dataclass(frozen=True)
class TileHeader:
    """What every tile of a run carries besides its units."""

    identity: str  # sha256 of the run's cfg.json
    grid: NDArray[np.float64]
    gap_every: int  # lambda_1 subsample of the units (0 = none)
    argv: tuple[str, ...]  # command line, for the provenance


class TileError(RuntimeError):
    """A tile file on disk that this run cannot use."""


class TileMismatch(TileError):
    """A tile file on disk belongs to another run, grid, tile or subsample."""


class TileUnreadable(TileError):
    """A tile file on disk is not a readable npz (empty, truncated, ...)."""


#: Arrays of every tile npz (see ``write_tile``).
TILE_KEYS = frozenset(
    (
        "schema",
        "identity",
        "provenance",
        "grid",
        "tile",
        "gap_every",
        "complete",
        "records",
        "ps_ptr",
        "ps_size",
        "ps_count",
    )
)
#: What ``np.load`` raises on a damaged npz (a missing file is not one).
_DAMAGED = (EOFError, OSError, ValueError, zipfile.BadZipFile, zlib.error)


def tile_stem(ghash: str, tile: Tile) -> str:
    """File stem of ``tile`` on the grid with hash ``ghash``."""
    return f"g={ghash}_c={tile.c0}-{tile.c1}_r={tile.r0}-{tile.r1}"


def tile_paths(rundir: Path, ghash: str, tile: Tile) -> TilePaths:
    """Where ``tile`` and its checkpoint live in ``rundir``."""
    stem = str(Path(rundir) / TC_TILES_DIR / tile_stem(ghash, tile))
    return TilePaths(Path(stem + TC_TILE_EXT), Path(stem + TC_PART_EXT))


class TileBuffer:
    """Finished units of one tile: record rows and their P(S)."""

    def __init__(self) -> None:
        self.records: list[np.void] = []
        self.ps: list[PS] = []

    def __len__(self) -> int:
        return len(self.records)

    def add(self, record: np.ndarray, ps: PS) -> None:
        self.records.append(record[()])
        self.ps.append(ps)

    def done(self) -> set[tuple[int, int]]:
        """The ``(rep, ictrl)`` of every finished unit."""
        return {(int(r["rep"]), int(r["ictrl"])) for r in self.records}

    def arrays(self) -> dict[str, NDArray]:
        """Records sorted by (rep, ictrl), with their P(S) as CSR arrays."""
        records = np.array(self.records, dtype=RECORD_DTYPE)
        order = np.lexsort((records["ictrl"], records["rep"]))
        sizes = [self.ps[i][0] for i in order]
        counts = [self.ps[i][1] for i in order]
        ptr = np.zeros(len(order) + 1, dtype=np.int64)
        ptr[1:] = np.cumsum([s.size for s in sizes], dtype=np.int64)
        return {
            "records": records[order],
            "ps_ptr": ptr,
            "ps_size": concat_int64(sizes),
            "ps_count": concat_int64(counts),
        }

    @classmethod
    def from_arrays(cls, data: dict[str, NDArray]) -> TileBuffer:
        """Inverse of ``arrays``."""
        buf, ptr = cls(), data["ps_ptr"]
        for i, record in enumerate(data["records"]):
            lo, hi = ptr[i], ptr[i + 1]
            buf.records.append(record)
            buf.ps.append((data["ps_size"][lo:hi], data["ps_count"][lo:hi]))
        return buf


def concat_int64(parts: Sequence[NDArray]) -> NDArray[np.int64]:
    """``parts`` joined as one int64 array (empty if there are none)."""
    if not parts:
        return np.zeros(0, dtype=np.int64)
    return np.concatenate(parts).astype(np.int64)


def write_tile(
    path: Path,
    buf: TileBuffer,
    header: TileHeader,
    tile: Tile,
    *,
    complete: bool,
) -> Path:
    """Atomically write ``buf`` as a tile npz (see the module docstring)."""
    return write_npz(path, _tile_arrays(buf, header, tile, complete=complete))


def _tile_arrays(
    buf: TileBuffer, header: TileHeader, tile: Tile, *, complete: bool
) -> dict[str, NDArray]:
    """The arrays of the tile npz of ``buf`` (see ``TILE_KEYS``)."""
    return {
        "schema": np.int64(TC_TILE_SCHEMA),
        "identity": np.str_(header.identity),
        "provenance": np.str_(provenance(header.argv)),
        "grid": header.grid,
        "tile": tile.bounds(),
        "gap_every": np.int64(header.gap_every),
        "complete": np.bool_(complete),
        **buf.arrays(),
    }


def write_npz(path: Path, arrays: dict[str, NDArray]) -> Path:
    """Atomically write ``arrays`` as a compressed npz."""
    data = _npz_bytes(arrays)
    return atomic_write_bytes(path, lambda fh: fh.write(data))


def _npz_bytes(arrays: dict[str, NDArray]) -> bytes:
    """``arrays`` as the bytes of a compressed npz."""
    payload = io.BytesIO()
    np.savez_compressed(payload, **arrays)
    return payload.getvalue()


def read_tile(path: Path) -> dict[str, NDArray]:
    """Every array of a tile npz (pickles refused).

    Raises ``TileUnreadable``, naming the file, if it is not a readable npz
    (for instance empty or truncated after a node crash).
    """
    path = Path(path)
    try:
        npz = np.load(path, allow_pickle=False)
        if not isinstance(npz, NpzFile):
            raise ValueError("not an npz archive")
        with npz:
            return {name: npz[name] for name in npz.files}
    except FileNotFoundError:
        raise
    except _DAMAGED as err:
        raise TileUnreadable(f"{path.name}: unreadable ({err})") from err


def load_tile(path: Path, identity: str) -> dict[str, NDArray]:
    """A complete tile of the run ``identity``, checked against its own
    grid and bounds. Raises ``TileError`` naming the file."""
    data = read_tile(path)
    failed = _tile_problems(data, identity, complete=True)
    if failed:
        raise _mismatch(path, failed)
    return data


def check_tile(
    path: Path,
    data: dict[str, NDArray],
    header: TileHeader,
    tile: Tile,
    *,
    complete: bool,
) -> None:
    """Raise ``TileMismatch``, naming ``path``, unless ``data`` is ``tile``
    of this run: same identity, grid, bounds and ``gap_every``. A complete
    tile must hold every unit of ``tile`` once, a partial a subset of them.
    """
    failed = _tile_problems(data, header.identity, complete=complete)
    if TILE_KEYS <= data.keys():
        run = {
            "grid": np.array_equal(data["grid"], header.grid),
            "bounds": np.array_equal(data["tile"], tile.bounds()),
            "gap_every": int(data["gap_every"]) == header.gap_every,
        }
        failed += [name for name, ok in run.items() if not ok]
    if failed:
        raise _mismatch(path, failed)


def _tile_problems(
    data: dict[str, NDArray], identity: str, *, complete: bool
) -> list[str]:
    """The checks that ``data`` fails as a tile of the run ``identity``,
    judged against the grid and bounds it stores."""
    missing = TILE_KEYS - data.keys()
    if missing:
        return [f"missing {', '.join(sorted(missing))}"]
    records, ptr = data["records"], data["ps_ptr"]
    checks = {
        "schema": int(data["schema"]) == TC_TILE_SCHEMA,
        "record dtype": records.dtype == RECORD_DTYPE,
        "identity": str(data["identity"]) == identity,
        "complete flag": bool(data["complete"]) == complete,
        "P(S) layout": ptr.size == records.size + 1
        and ptr[0] == 0
        and bool(np.all(np.diff(ptr) >= 0))
        and ptr[-1] == data["ps_size"].size == data["ps_count"].size,
    }
    if checks["record dtype"]:
        tile = Tile(*data["tile"].tolist())
        keys = set(zip(records["rep"].tolist(), records["ictrl"].tolist()))
        wanted = tile.units()
        covered = keys == wanted if complete else keys <= wanted
        checks["units"] = covered and len(keys) == len(records)
    return [name for name, ok in checks.items() if not ok]


def _mismatch(path: Path, failed: Sequence[str]) -> TileMismatch:
    names = ", ".join(failed)
    return TileMismatch(f"{Path(path).name}: not a tile of this run ({names})")


def set_aside(path: Path, reason: Exception) -> Optional[Path]:
    """Rename a bad tile file to ``<name>.<id>.bad`` and warn; return the
    new path (``None`` if another task already moved it)."""
    token = secrets.token_hex(TC_BAD_ID_BYTES)
    bad = path.with_name(f"{path.name}.{token}{TC_BAD_EXT}")
    try:
        os.replace(path, bad)
    except FileNotFoundError:
        return None
    warnings.warn(f"{reason}; renamed to {bad.name}", stacklevel=2)
    return bad


def resume_buffer(part: Path, header: TileHeader, tile: Tile) -> TileBuffer:
    """Units already finished in the checkpoint of ``tile`` (none if absent).

    A checkpoint that is unreadable or belongs to another run is set aside,
    and the tile starts over.
    """
    if not part.exists():
        return TileBuffer()
    try:
        data = read_tile(part)
        check_tile(part, data, header, tile, complete=False)
    except TileError as err:
        set_aside(part, err)
        return TileBuffer()
    return TileBuffer.from_arrays(data)


def _final_done(final: Path, header: TileHeader, tile: Tile) -> bool:
    """Whether ``final`` already holds ``tile``. An unreadable final is set
    aside (to be recomputed); a final of another run raises."""
    if not final.exists():
        return False
    try:
        data = read_tile(final)
    except TileUnreadable as err:
        set_aside(final, err)
        return False
    check_tile(final, data, header, tile, complete=True)
    return True


def publish_tile(
    final: Path, buf: TileBuffer, header: TileHeader, tile: Tile
) -> None:
    """Hard-link ``buf`` onto ``final`` as the complete ``tile``.

    A file at ``final`` is never replaced. One published meanwhile by another
    task is kept if it holds ``tile`` of this run and raises ``TileMismatch``
    otherwise; an unreadable one is set aside and the link retried once.
    """
    data = _npz_bytes(_tile_arrays(buf, header, tile, complete=True))
    if create_once(final, data) or _final_done(final, header, tile):
        return
    if not create_once(final, data):  # another final since the set-aside
        check_tile(final, read_tile(final), header, tile, complete=True)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
def provenance(argv: Sequence[str]) -> str:
    """JSON provenance of this process: command, SLURM ids, host, versions.

    Taken at write time, so it lists the modules the units imported.
    """
    modules = {name: _module_version(name) for name in TC_PROVENANCE_MODULES}
    info = {
        "argv": shlex.join(argv),
        "slurm": {name: os.environ.get(name) for name in TC_SLURM_ENV_VARS},
        "host": socket.gethostname(),
        "cpu": _cpu_model(),
        "python": platform.python_version(),
        "modules": modules,
        "git": git_state(),
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    return json.dumps(info, sort_keys=True)


def _module_version(name: str) -> Optional[str]:
    """Version of an already imported module (never imports it)."""
    return getattr(sys.modules.get(name), "__version__", None)


def _cpu_model() -> str:
    """CPU model from ``/proc/cpuinfo``, else the platform's guess."""
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


@functools.lru_cache(maxsize=1)
def git_state() -> Optional[dict[str, Any]]:
    """Commit and dirty flag of the checkout that holds these programs.

    The first call is cached: programs call it at start, so it describes
    the code they imported. Untracked files under ``src/`` count as dirty.
    Git runs without optional locks and never writes the checkout's index.
    """
    root = Path(__file__).resolve().parents[1]
    try:
        commit = _git(root, "rev-parse", "HEAD")
        status = _git(
            root, "status", "--porcelain", "--untracked-files=normal", "--", "."
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return {"commit": commit, "dirty": bool(status)}


def _git(root: Path, *args: str) -> str:
    out = subprocess.run(
        ["git", "--no-optional-locks", "-C", str(root), *args],
        capture_output=True,
        text=True,
        check=True,
        timeout=TC_GIT_TIMEOUT,
    )
    return out.stdout.strip()


# ---------------------------------------------------------------------------
# The tile loop
# ---------------------------------------------------------------------------
class StopFlag:
    """Remembers the last stop signal received (``None`` if none)."""

    def __init__(self) -> None:
        self.signum: Optional[int] = None

    def __call__(self, signum: int, frame: object) -> None:
        self.signum = signum


@contextmanager
def stop_signals() -> Iterator[StopFlag]:
    """Route ``TC_STOP_SIGNALS`` to a flag; restore the handlers on exit."""
    flag = StopFlag()
    signums = [getattr(signal, name, None) for name in TC_STOP_SIGNALS]
    previous = {s: signal.signal(s, flag) for s in signums if s is not None}
    try:
        yield flag
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def pending_units(
    tile: Tile, done: set[tuple[int, int]]
) -> Iterator[tuple[int, list[int]]]:
    """``(rep, ictrls)`` still to run in ``tile``, in ascending rep order."""
    for rep in range(tile.r0, tile.r1):
        ictrls = [c for c in range(tile.c0, tile.c1) if (rep, c) not in done]
        if ictrls:
            yield rep, ictrls


def run_tile(
    sweep: Sweep,
    tile: Tile,
    paths: TilePaths,
    header: TileHeader,
    stop: StopFlag,
    ckpt_seconds: float = TC_CKPT_SECONDS,
) -> int:
    """Run the missing units of ``tile`` and return an exit code.

    A complete tile on disk makes this a no-op, and a checkpoint is resumed
    (bad files are set aside first, see the module docstring). Once ``stop``
    holds a signal, the next unit boundary writes the checkpoint and returns
    ``TC_EXIT_PREEMPTED``; the last unit publishes the complete tile instead
    (``publish_tile``), and the caller checks ``stop`` before the next tile.
    Raises ``TileMismatch`` if a complete tile of another run is on disk, at
    the start or when the tile is published.
    """
    if _final_done(paths.final, header, tile):
        return TC_EXIT_OK
    buf = resume_buffer(paths.part, header, tile)
    todo = pending_units(tile, buf.done())
    units = itertools.chain.from_iterable(sweep(r, c) for r, c in todo)
    last = time.monotonic()
    for record, ps in units:
        buf.add(record, ps)
        stopping = stop.signum is not None
        due = time.monotonic() - last >= ckpt_seconds
        if len(buf) == tile.size or not (stopping or due):
            continue
        write_tile(paths.part, buf, header, tile, complete=False)
        last = time.monotonic()
        if stopping:
            return TC_EXIT_PREEMPTED
    if buf.done() != tile.units():
        raise RuntimeError("the sweep did not return the tile's units")
    publish_tile(paths.final, buf, header, tile)
    paths.part.unlink(missing_ok=True)
    return TC_EXIT_OK
