"""RRG_TransCluster program: tiling, resume, run identity and command line."""

import dataclasses
import json
import os
import random
import re
import signal
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[2] / "src"
sys.path.insert(0, str(SRC))

import kernels.RRG_TransCluster as glue  # noqa: E402
import kernels.TransClusterTiles as tiles  # noqa: E402
from kernels.TransClusterCore import RECORD_DTYPE, make_record  # noqa: E402
from kernels.TransClusterTiles import (  # noqa: E402
    TILE_KEYS,
    StopFlag,
    Tile,
    TileBuffer,
    TileHeader,
    all_tiles,
    decode_tile,
    grid_hash,
    make_layout,
    read_tile,
    resolve_grid,
    run_tile,
    tile_paths,
    write_tile,
)
from lrgsglib.config.progargs.defs.TransCluster import (  # noqa: E402
    DEFAULT_TC_GAP_EVERY,
    TC_BAD_EXT,
    TC_BAD_ID_BYTES,
    TC_EIG_TOL_MAX,
    TC_EXIT_CFG_MISMATCH,
    TC_EXIT_OK,
    TC_EXIT_PREEMPTED,
    TC_PART_EXT,
    TC_RUNDIR_NAME_MAX,
    TC_SEED_SPAN,
    TC_STOP_SIGNALS,
    TC_TILES_DIR,
)
from lrgsglib.statsys._io import publish_once  # noqa: E402
from lrgsglib.statsys._naming import cfg_sidecar_bytes  # noqa: E402
from parsers.RRG_TransCluster import build_parser  # noqa: E402
from RRG_TransCluster import main  # noqa: E402

pytestmark = pytest.mark.integration

N, K, REPS, N_CTRL = 60, 3, 6, 6
GRID = f"(0, 0.5, {N_CTRL})"
BASE = [str(N), str(K), "--p-linsp", GRID, "-na", str(REPS)]
N_UNITS = N_CTRL * REPS
REAL_SWEEP = glue.sweep_rep
TILINGS = [("0", "0"), ("1", "1"), ("4", "5")]
PUBLISHERS = 16
STOP_SIGNUMS = [getattr(signal, name) for name in TC_STOP_SIGNALS]
# A 3-control x 2-rep tile of synthetic units, for the bare tile loop
SYN_TILE = decode_tile(make_layout(3, 0, 2, 0, 0), 0)
SYN_HEADER = TileHeader("id", np.linspace(0.0, 0.2, 3), 0, ("prog",))
# -os tags whose run directory name has more than TC_RUNDIR_NAME_MAX bytes;
# the wide one has fewer characters than that (U+00E9 takes two bytes)
LONG_TAG = "a" * TC_RUNDIR_NAME_MAX
WIDE_TAG = "\u00e9" * (TC_RUNDIR_NAME_MAX // 2)
ARGUMENT_ERRORS = {
    "k below the minimum": "60 2 --p-list 0.1",
    "N*k odd": "61 3 --p-list 0.1",
    "no grid": "60 3",
    "two grids": "60 3 --p-list 0.1 --sigma-list 1",
    "sigma grid for flip": "60 3 --sigma-list 0.5",
    "p grid for gauss": "60 3 -ew normal --p-list 0.1",
    "mu for flip": "60 3 --mu 0.5 --p-list 0.1",
    "ensemble for gauss": "60 3 -ew normal --ensemble bernoulli --sigma-list 1",
    "cell for gauss": "60 3 -ew normal -c randXERR --sigma-list 1",
    "p above one": "60 3 --p-list 1.5",
    "sigma 0 at mu 0": "60 3 -ew normal --mu 0 --sigma-list 0 0.5",
    "coupling scale below the window": (
        "60 3 -ew normal --mu 0 --sigma-list 1e-26"
    ),
    "sigma above the scale window": "60 3 -ew normal --sigma-list 1e101",
    "mu below the scale window": "60 3 -ew normal --mu 1e-9 --sigma-list 0",
    "mu above the scale window": "60 3 -ew normal --mu 1e200 --sigma-list 0.5",
    "NaN mu": "60 3 -ew normal --mu nan --sigma-list 0.5",
    "infinite mu": "60 3 -ew normal --mu inf --sigma-list 0.5",
    "tile id out of range": "60 3 --p-list 0.1 -na 6 --rep-chunk 0 -tid 1 -ge nx",
    "tag with separator": "60 3 --p-list 0.1 -os a_b",
    "run dir name too long": f"60 3 --p-list 0.1 -os {LONG_TAG}",
    "run dir name too long in bytes": f"60 3 --p-list 0.1 -os {WIDE_TAG}",
    "negative seed": "60 3 --p-list 0.1 --seed -1",
    "reps past the seed span": (
        f"60 3 --p-list 0.1 --rep-start {TC_SEED_SPAN} -na 1"
    ),
    "NaN eig-tol": "60 3 --p-list 0.1 --eig-tol nan",
    "infinite eig-tol": "60 3 --p-list 0.1 --eig-tol inf",
    "negative eig-tol": "60 3 --p-list 0.1 --eig-tol -1",
    "eig-tol typo 1e6": "60 3 --p-list 0.1 --eig-tol 1e6",
    "eig-tol above the cap": (
        f"60 3 --p-list 0.1 --eig-tol {2 * TC_EIG_TOL_MAX}"
    ),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _run(root: Path, monkeypatch, *extra: str) -> int:
    pytest.importorskip("graph_tool")  # the program's default engine
    monkeypatch.setattr(glue, "PATHDATA", root)
    return main([*BASE, *extra])


def _resolve(argv: list[str]) -> glue.RRGRun:
    """The run of ``argv`` on the default engine (needs graph-tool)."""
    pytest.importorskip("graph_tool")
    return glue.resolve_run(build_parser().parse_args(argv), argv)


def _rundir(root: Path) -> Path:
    (rundir,) = root.glob(f"*/phtra/N={N}_k={K}/*")
    return rundir


def _tile_files(root: Path, pattern: str = "*") -> list[Path]:
    return sorted((_rundir(root) / TC_TILES_DIR).glob(pattern))


def _records(rundir: Path) -> tuple[np.ndarray, list]:
    """Every record of the complete tiles, sorted by (rep, ictrl)."""
    records, ps = [], []
    for path in sorted((rundir / TC_TILES_DIR).glob("*.npz")):
        if not path.name.endswith(TC_PART_EXT):
            buf = TileBuffer.from_arrays(read_tile(path))
            records, ps = records + buf.records, ps + buf.ps
    table = np.array(records, dtype=RECORD_DTYPE)
    order = np.lexsort((table["ictrl"], table["rep"]))
    return table[order], [ps[i] for i in order]


def _assert_table_equal(a: np.ndarray, b: np.ndarray, skip=()) -> None:
    """Field-wise equality (``==`` on structured rows treats NaN as unequal)."""
    for name in a.dtype.names:
        if name not in skip:
            np.testing.assert_array_equal(a[name], b[name], name)


def _assert_same_records(a: tuple, b: tuple) -> None:
    """Bitwise equal records and P(S), except the solve time."""
    (rec_a, ps_a), (rec_b, ps_b) = a, b
    assert len(rec_a) == len(rec_b) == N_UNITS
    _assert_table_equal(rec_a, rec_b, skip=("t_eig",))
    for (size_a, count_a), (size_b, count_b) in zip(ps_a, ps_b):
        np.testing.assert_array_equal(size_a, size_b)
        np.testing.assert_array_equal(count_a, count_b)


def _counting(monkeypatch, *, stop_after=0, signum=signal.SIGUSR1) -> list:
    """Wrap the core sweep; count units and optionally signal ourselves."""
    produced = []

    def sweep(*args, **kwargs):
        for unit in REAL_SWEEP(*args, **kwargs):
            produced.append(unit)
            if len(produced) == stop_after:
                signal.raise_signal(signum)  # handled before the yield
            yield unit

    monkeypatch.setattr(glue, "sweep_rep", sweep)
    return produced


def _published_meanwhile(monkeypatch, source: Path, root: Path) -> None:
    """Wrap the core sweep: once the units start (after run_tile's start
    check), another task publishes ``source`` as the complete tile."""

    def sweep(*args, **kwargs):
        final = _rundir(root) / TC_TILES_DIR / source.name
        if not final.exists():
            final.parent.mkdir(exist_ok=True)
            final.write_bytes(source.read_bytes())
        yield from REAL_SWEEP(*args, **kwargs)

    monkeypatch.setattr(glue, "sweep_rep", sweep)


def _handlers() -> dict:
    return {signum: signal.getsignal(signum) for signum in STOP_SIGNUMS}


def _rewrite_header(path: Path, **changes) -> None:
    """Rewrite a tile file with some of its header fields changed."""
    data = read_tile(path)
    header = TileHeader(
        str(data["identity"]), data["grid"], int(data["gap_every"]), ("prog",)
    )
    header = dataclasses.replace(header, **changes)
    tile, complete = Tile(*data["tile"].tolist()), bool(data["complete"])
    write_tile(
        path, TileBuffer.from_arrays(data), header, tile, complete=complete
    )


def _truncate(path: Path) -> None:
    path.write_bytes(path.read_bytes()[: path.stat().st_size // 2])


def _synthetic_sweep(rep, ictrls):
    for ictrl in ictrls:
        yield make_record({"rep": rep, "ictrl": ictrl}), (
            np.array([1, 2], dtype=np.int64),
            np.array([rep + 1, ictrl + 1], dtype=np.int64),
        )


@pytest.fixture
def stop_guard():
    """Harmless stop-signal handlers, so a missed handler cannot kill pytest."""
    previous = {s: signal.signal(s, lambda *_: None) for s in STOP_SIGNUMS}
    yield
    for signum, handler in previous.items():
        signal.signal(signum, handler)


# ---------------------------------------------------------------------------
# Grid and layout
# ---------------------------------------------------------------------------
def test_grid_is_canonical_sorted_and_hashed():
    grid = resolve_grid([0.3, 0.1 + 0.2, 0.1], 0.0, 1.0)
    np.testing.assert_array_equal(grid, [0.1, 0.3])
    assert grid_hash(grid) == grid_hash(resolve_grid([0.1, 0.3], 0.0, 1.0))
    assert grid_hash(grid) != grid_hash(resolve_grid([0.1, 0.31], 0.0, 1.0))
    zero = resolve_grid([-0.0, 0.1], 0.0, 1.0)  # -0.0 is the point 0.0
    assert grid_hash(zero) == grid_hash(resolve_grid([0.0, 0.1], 0.0, 1.0))
    for bad in ([], [0.5, 1.5], [-0.1], [float("nan")]):
        with pytest.raises(ValueError):
            resolve_grid(bad, 0.0, 1.0)


@pytest.mark.parametrize(
    "n_ctrl,n_reps,ctrl_chunk,rep_chunk,rep_start",
    [
        (7, 6, 3, 4, 0),
        (7, 6, 0, 0, 0),
        (7, 6, 1, 1, 5),
        (26, 200, 5, 7, 1000),
        (1, 1, 0, 0, 0),
        (5, 3, 9, 9, 0),
    ],
)
def test_tiles_cover_every_unit_once(
    n_ctrl, n_reps, ctrl_chunk, rep_chunk, rep_start
):
    layout = make_layout(n_ctrl, rep_start, n_reps, ctrl_chunk, rep_chunk)
    units = [u for tile in all_tiles(layout) for u in tile.units()]
    reps = range(rep_start, rep_start + n_reps)
    assert len(units) == len(set(units)) == n_ctrl * n_reps
    assert set(units) == set(product(reps, range(n_ctrl)))


def test_tile_ids_decode_rep_major():
    layout = make_layout(7, 0, 6, 3, 4)
    first = [decode_tile(layout, t) for t in range(layout.n_ctrl_tiles)]
    assert {(t.r0, t.r1) for t in first} == {(0, 4)}
    assert [(t.c0, t.c1) for t in first] == [(0, 3), (3, 6), (6, 7)]
    assert decode_tile(layout, layout.n_ctrl_tiles).r0 == 4
    for bad in (-1, layout.n_tiles):
        with pytest.raises(ValueError):
            decode_tile(layout, bad)


def test_tile_buffer_round_trip():
    buf = TileBuffer()
    for rep, ictrl in [(1, 0), (0, 1), (0, 0)]:
        for rec, ps in _synthetic_sweep(rep, [ictrl]):
            buf.add(rec, ps)
    arrays = buf.arrays()
    keys = arrays["records"][["rep", "ictrl"]].tolist()
    assert keys == [(0, 0), (0, 1), (1, 0)]
    again = TileBuffer.from_arrays(arrays).arrays()
    assert arrays.keys() == again.keys()
    _assert_table_equal(arrays.pop("records"), again.pop("records"))
    for name, value in arrays.items():
        np.testing.assert_array_equal(value, again[name], name)


# ---------------------------------------------------------------------------
# The tile loop
# ---------------------------------------------------------------------------
def test_checkpoints_follow_the_clock(tmp_path, monkeypatch):
    tile, header = SYN_TILE, SYN_HEADER
    writes = []
    real_write, real_publish = tiles.write_tile, tiles.publish_tile

    def spy(path, buf, header, tile, *, complete):
        writes.append((len(buf), complete))
        return real_write(path, buf, header, tile, complete=complete)

    def spy_publish(final, buf, header, tile):
        writes.append((len(buf), True))
        return real_publish(final, buf, header, tile)

    monkeypatch.setattr(tiles, "write_tile", spy)
    monkeypatch.setattr(tiles, "publish_tile", spy_publish)
    paths = tile_paths(tmp_path, "g", tile)
    stop = StopFlag()
    code = run_tile(_synthetic_sweep, tile, paths, header, stop, ckpt_seconds=0)
    assert code == TC_EXIT_OK
    assert writes == [(n, False) for n in range(1, 6)] + [(6, True)]
    assert paths.final.exists() and not paths.part.exists()
    assert read_tile(paths.final).keys() == TILE_KEYS


def test_a_raised_stop_flag_checkpoints_at_the_next_unit(tmp_path):
    tile, header = SYN_TILE, SYN_HEADER
    paths, stop = tile_paths(tmp_path, "g", tile), StopFlag()
    stop(signal.SIGTERM, None)
    code = run_tile(_synthetic_sweep, tile, paths, header, stop)
    assert code == TC_EXIT_PREEMPTED and not paths.final.exists()
    part = read_tile(paths.part)
    assert part.keys() == TILE_KEYS and len(part["records"]) == 1
    code = run_tile(_synthetic_sweep, tile, paths, header, StopFlag())
    assert code == TC_EXIT_OK and not paths.part.exists()
    assert len(read_tile(paths.final)["records"]) == tile.size


def test_records_do_not_depend_on_tiling_or_order(tmp_path, monkeypatch):
    runs = []
    for ctrl_chunk, rep_chunk in TILINGS:
        root = tmp_path / f"c{ctrl_chunk}_r{rep_chunk}"
        chunks = ["--ctrl-chunk", ctrl_chunk, "--rep-chunk", rep_chunk]
        assert _run(root, monkeypatch, *chunks) == TC_EXIT_OK
        runs.append(_records(_rundir(root)))
    layout = make_layout(N_CTRL, 0, REPS, 4, 5)
    order = list(range(layout.n_tiles))
    random.Random(0).shuffle(order)
    root = tmp_path / "shuffled"
    for tile_id in order:
        extra = ["--ctrl-chunk", "4", "--rep-chunk", "5", "-tid", str(tile_id)]
        assert _run(root, monkeypatch, *extra) == TC_EXIT_OK
    runs.append(_records(_rundir(root)))
    assert any(ps[0].size for ps in runs[0][1]), "P(S) never exercised"
    for other in runs[1:]:
        _assert_same_records(runs[0], other)


def test_complete_tiles_are_not_recomputed(tmp_path, monkeypatch):
    assert _run(tmp_path, monkeypatch) == TC_EXIT_OK
    produced = _counting(monkeypatch)
    assert _run(tmp_path, monkeypatch) == TC_EXIT_OK
    assert produced == []


@pytest.mark.parametrize("signame", TC_STOP_SIGNALS)
def test_stop_signal_checkpoints_and_resumes(
    signame, tmp_path, monkeypatch, stop_guard
):
    handlers = _handlers()
    ref = tmp_path / "ref"
    assert _run(ref, monkeypatch, "--rep-chunk", "0") == TC_EXIT_OK
    reference = _records(_rundir(ref))
    root = tmp_path / "cut"
    _counting(monkeypatch, stop_after=3, signum=getattr(signal, signame))
    code = _run(root, monkeypatch, "--rep-chunk", "0")
    assert code == TC_EXIT_PREEMPTED
    assert _handlers() == handlers
    (part,) = _tile_files(root)
    assert part.name.endswith(TC_PART_EXT)
    assert len(read_tile(part)["records"]) == 3
    produced = _counting(monkeypatch)
    assert _run(root, monkeypatch, "--rep-chunk", "0") == TC_EXIT_OK
    assert len(produced) == N_UNITS - 3
    assert not part.exists() and _handlers() == handlers
    _assert_same_records(reference, _records(_rundir(root)))


def test_a_stop_in_the_last_unit_of_a_tile_ends_the_run_after_it(
    tmp_path, monkeypatch, stop_guard
):
    """-tid -1: the signalled tile completes, the next one never starts."""
    chunks, per_tile = ("--rep-chunk", str(REPS // 2)), N_UNITS // 2
    ref = tmp_path / "ref"
    assert _run(ref, monkeypatch, *chunks) == TC_EXIT_OK
    reference = _records(_rundir(ref))
    root = tmp_path / "cut"
    produced = _counting(monkeypatch, stop_after=per_tile)
    assert _run(root, monkeypatch, *chunks) == TC_EXIT_PREEMPTED
    assert len(produced) == per_tile
    (first,) = _tile_files(root)
    assert not first.name.endswith(TC_PART_EXT)
    produced = _counting(monkeypatch)
    assert _run(root, monkeypatch, *chunks) == TC_EXIT_OK
    assert len(produced) == per_tile
    _assert_same_records(reference, _records(_rundir(root)))


BAD_TILES = {  # damage, and the units run before it (0: a complete tile)
    "empty final": (lambda path: path.write_bytes(b""), 0),
    "truncated checkpoint": (_truncate, 3),
    "checkpoint of another gap_every": (
        lambda path: _rewrite_header(path, gap_every=DEFAULT_TC_GAP_EVERY + 1),
        3,
    ),
}


@pytest.mark.parametrize("damage,stop_after", BAD_TILES.values(), ids=BAD_TILES)
def test_bad_tile_files_are_set_aside_and_recomputed(
    damage, stop_after, tmp_path, monkeypatch, stop_guard
):
    ref = tmp_path / "ref"
    assert _run(ref, monkeypatch, "--rep-chunk", "0") == TC_EXIT_OK
    reference = _records(_rundir(ref))
    root = tmp_path / "bad"
    _counting(monkeypatch, stop_after=stop_after)
    _run(root, monkeypatch, "--rep-chunk", "0")
    (path,) = _tile_files(root)
    damage(path)
    produced = _counting(monkeypatch)
    with pytest.warns(UserWarning, match="renamed to"):
        assert _run(root, monkeypatch, "--rep-chunk", "0") == TC_EXIT_OK
    assert len(produced) == N_UNITS  # nothing reused from the bad file
    (bad,) = _tile_files(root, f"*{TC_BAD_EXT}")
    token = f"[0-9a-f]{{{2 * TC_BAD_ID_BYTES}}}"
    name, ext = re.escape(path.name), re.escape(TC_BAD_EXT)
    assert re.fullmatch(rf"{name}\.{token}{ext}", bad.name)
    (final,) = _tile_files(root, "*.npz")
    assert not final.name.endswith(TC_PART_EXT)
    _assert_same_records(reference, _records(_rundir(root)))


RACES = {  # --gap-every of the task that publishes first, and our exit code
    "same run": (DEFAULT_TC_GAP_EVERY, TC_EXIT_OK),
    "another gap_every": (DEFAULT_TC_GAP_EVERY + 1, TC_EXIT_CFG_MISMATCH),
}


@pytest.mark.parametrize("gap_every,code", RACES.values(), ids=RACES)
def test_a_final_published_meanwhile_is_never_replaced(
    gap_every, code, tmp_path, monkeypatch, stop_guard
):
    """Another task completes the tile after our start check. Its final is
    kept; our checkpoint goes if it is a tile of our run, else we exit 3."""
    chunks = ("--rep-chunk", "0")
    other = tmp_path / "other"
    gap = ("--gap-every", str(gap_every))
    assert _run(other, monkeypatch, *chunks, *gap) == TC_EXIT_OK
    (theirs,) = _tile_files(other)
    root = tmp_path / "race"
    _counting(monkeypatch, stop_after=3)
    assert _run(root, monkeypatch, *chunks) == TC_EXIT_PREEMPTED
    (part,) = _tile_files(root)
    _published_meanwhile(monkeypatch, theirs, root)
    assert _run(root, monkeypatch, *chunks) == code
    final = part.parent / theirs.name
    assert final.read_bytes() == theirs.read_bytes()
    left = [final] if code == TC_EXIT_OK else [final, part]
    assert _tile_files(root) == sorted(left)  # and no temp file


def test_an_unreadable_final_published_meanwhile_is_set_aside(
    tmp_path, monkeypatch
):
    ref = tmp_path / "ref"
    assert _run(ref, monkeypatch, "--rep-chunk", "0") == TC_EXIT_OK
    (ref_tile,) = _tile_files(ref)
    empty = tmp_path / ref_tile.name
    empty.write_bytes(b"")
    root = tmp_path / "race"
    _published_meanwhile(monkeypatch, empty, root)
    with pytest.warns(UserWarning, match="renamed to"):
        assert _run(root, monkeypatch, "--rep-chunk", "0") == TC_EXIT_OK
    (bad,) = _tile_files(root, f"*{TC_BAD_EXT}")
    assert bad.read_bytes() == b"" and len(_tile_files(root)) == 2
    _assert_same_records(_records(_rundir(ref)), _records(_rundir(root)))


# ---------------------------------------------------------------------------
# Run identity
# ---------------------------------------------------------------------------
def test_foreign_cfg_exits_with_mismatch(tmp_path, monkeypatch, capsys):
    assert _run(tmp_path, monkeypatch) == TC_EXIT_OK
    cfg = _rundir(tmp_path) / "cfg.json"
    payload = json.loads(cfg.read_text())
    payload["graph"]["seed"] += 1
    cfg.write_text(json.dumps(payload))
    handlers = _handlers()
    assert _run(tmp_path, monkeypatch) == TC_EXIT_CFG_MISMATCH
    assert "different content" in capsys.readouterr().err
    assert _handlers() == handlers


def test_foreign_tile_exits_with_mismatch(tmp_path, monkeypatch, capsys):
    assert _run(tmp_path, monkeypatch, "--rep-chunk", "0") == TC_EXIT_OK
    (path,) = _tile_files(tmp_path)
    _rewrite_header(path, identity="another run")
    code = _run(tmp_path, monkeypatch, "--rep-chunk", "0")
    assert code == TC_EXIT_CFG_MISMATCH
    assert "identity" in capsys.readouterr().err
    assert _tile_files(tmp_path) == [path]  # an error: left in place


def test_complete_tile_of_another_gap_every_exits_with_mismatch(
    tmp_path, monkeypatch, capsys
):
    assert _run(tmp_path, monkeypatch) == TC_EXIT_OK
    produced = _counting(monkeypatch)
    other = str(DEFAULT_TC_GAP_EVERY + 1)
    code = _run(tmp_path, monkeypatch, "--gap-every", other)
    assert code == TC_EXIT_CFG_MISMATCH
    assert "gap_every" in capsys.readouterr().err
    assert produced == []


@pytest.mark.parametrize("flag", ["--eig-tol", "--mu"])
def test_a_negative_zero_is_the_run_of_zero(flag, tmp_path, monkeypatch):
    """-0.0 names the run directory and cfg.json of 0.0: the tile is reused."""
    pytest.importorskip("graph_tool")
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    argv = [str(N), str(K), "-ew", "normal", "--sigma-list", "0.5", "-na", "2"]
    assert main([*argv, flag, "0"]) == TC_EXIT_OK
    produced = _counting(monkeypatch)
    assert main([*argv, flag, "-0.0"]) == TC_EXIT_OK
    assert produced == []
    _rundir(tmp_path)  # a single run directory


def test_reps_may_reach_the_seed_span(tmp_path, monkeypatch):
    last = TC_SEED_SPAN - 1
    extra = ["-na", "1", "--rep-start", str(last)]
    assert _run(tmp_path, monkeypatch, *extra) == TC_EXIT_OK
    records, _ = _records(_rundir(tmp_path))
    assert records["rep"].tolist() == [last] * N_CTRL


def test_publish_once_has_one_winner(tmp_path):
    """Racing publishers: os.link lets exactly one of them create the file."""
    path, data = tmp_path / "cfg.json", b'{"run": 1}\n'
    start = threading.Barrier(PUBLISHERS)

    def publish(_):
        start.wait()
        return publish_once(path, data)

    with ThreadPoolExecutor(PUBLISHERS) as pool:
        created = list(pool.map(publish, range(PUBLISHERS)))
    assert sum(created) == 1
    assert path.read_bytes() == data
    assert not [p for p in tmp_path.iterdir() if p != path]
    with pytest.raises(FileExistsError):
        publish_once(path, b'{"run": 2}\n')


@pytest.mark.parametrize("value", [float("nan"), float("inf")])
def test_cfg_json_is_strict_json(value):
    with pytest.raises(ValueError):
        cfg_sidecar_bytes(graph={}, model={"mu": value}, run={})


def test_a_non_json_identity_exits_through_argparse(tmp_path, monkeypatch):
    """A NaN in the identity is an argument error, not a crash at publish."""
    pytest.importorskip("graph_tool")  # the program's default engine
    real = glue.identity_blocks

    def nan_identity(*args):
        blocks = real(*args)
        return {**blocks, "model": {**blocks["model"], "mu": float("nan")}}

    monkeypatch.setattr(glue, "identity_blocks", nan_identity)
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    with pytest.raises(SystemExit) as exc:
        main(BASE)
    assert exc.value.code == 2
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize(
    "extra,name",
    [
        ([], "lang=np_s=20260930"),
        (
            ["--eig-tol", "1e-8", "-c", "randXERR", "--ensemble", "bernoulli"],
            "tol=1e-08_cell=randXERR_ens=bernoulli_lang=np_s=20260930",
        ),
        (["--seed", "7", "-os", "pilot"], "tag=pilot_lang=np_s=7"),
        (["--sampler", "native"], "smp=native_lang=np_s=20260930"),
    ],
)
def test_run_directory_names_the_identity(extra, name):
    run = _resolve([*BASE, *extra])
    assert run.rundir.name == name
    assert run.rundir.parent.name == f"N={N}_k={K}"


def test_gaussian_run_directory():
    run = _resolve(
        ["60", "3", "-ew", "normal", "--mu", "0.5", "--sigma-list", "1"]
    )
    assert run.rundir.name == "mu=0.5_ew=gauss_lang=np_s=20260930"
    assert run.identity["model"] == {"family": "G", "ensemble": None, "mu": 0.5}


def test_native_nx_warns_and_names_its_engine():
    argv = [*BASE, "--sampler", "native", "-ge", "nx"]
    with pytest.warns(UserWarning, match="not uniform"):
        run = glue.resolve_run(build_parser().parse_args(argv), argv)
    assert run.rundir.name == "smp=native_ge=nx_lang=np_s=20260930"
    assert run.identity["graph"]["engine"] == "nx"


@pytest.mark.parametrize("argv", ARGUMENT_ERRORS.values(), ids=ARGUMENT_ERRORS)
def test_argument_errors_exit_through_argparse(argv, tmp_path, monkeypatch):
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    with pytest.raises(SystemExit) as exc:
        main(argv.split())
    assert exc.value.code == 2
    assert not any(tmp_path.iterdir())


# ---------------------------------------------------------------------------
# Engines
# ---------------------------------------------------------------------------
def test_pairing_gives_the_same_records_on_both_engines(tmp_path, monkeypatch):
    runs = []
    for engine in ("nx", "gt"):
        root = tmp_path / engine
        assert _run(root, monkeypatch, "-ge", engine) == TC_EXIT_OK
        runs.append(_records(_rundir(root)))
    _assert_same_records(*runs)
    assert _rundir(tmp_path / "nx").name == _rundir(tmp_path / "gt").name
    produced = _counting(monkeypatch)  # one ensemble: gt reuses nx's tiles
    assert _run(tmp_path / "nx", monkeypatch, "-ge", "gt") == TC_EXIT_OK
    assert produced == []


@pytest.mark.slow
def test_native_gt_is_deterministic_across_processes(tmp_path):
    pytest.importorskip("graph_tool")
    runs = []
    for label in ("a", "b"):
        env = {**os.environ, "LRGSG_DATA": str(tmp_path / label)}
        command = [sys.executable, str(SRC / "RRG_TransCluster.py"), *BASE]
        command += ["--sampler", "native", "-ge", "gt"]
        done = subprocess.run(command, env=env, capture_output=True, text=True)
        assert done.returncode == TC_EXIT_OK, done.stderr
        runs.append(_records(_rundir(tmp_path / label)))
    _assert_same_records(*runs)
    first = next((_rundir(tmp_path / "a") / TC_TILES_DIR).glob("*"))
    provenance = json.loads(str(read_tile(first)["provenance"]))
    assert provenance["modules"]["graph_tool"]
    assert "--sampler native" in provenance["argv"]
