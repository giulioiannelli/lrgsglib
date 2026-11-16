"""RRG_TransCluster_Serializer: arrays, status, resubmission and collect."""

import dataclasses
import errno
import json
import os
import random
import re
import shlex
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("graph_tool")  # the program's default engine

SRC = Path(__file__).resolve().parents[2] / "src"
sys.path.insert(0, str(SRC))

import kernels.RRG_TransCluster as glue  # noqa: E402
import kernels.RRG_TransCluster_Serializer as ser  # noqa: E402
from kernels.TransClusterStats import load_agg  # noqa: E402
from kernels.TransClusterTiles import (  # noqa: E402
    TilePaths,
    decode_tile,
    grid_hash,
    read_tile,
    tile_paths,
    write_npz,
)
from lrgsglib.config.progargs.defs.generic import SBATCH_LOG_DIR  # noqa: E402
from lrgsglib.config.progargs.defs.TransCluster import (  # noqa: E402
    TC_AGG_NAME,
    TC_BAD_EXT,
    TC_BOOT_MIN_REPS,
    TC_EXIT_CFG_MISMATCH,
    TC_EXIT_FAILED,
    TC_EXIT_OK,
    TC_SUBMISSIONS_DIR,
)
from lrgsglib.statsys._naming import (  # noqa: E402
    CFG_SIDECAR_NAME,
    publish_cfg_sidecar,
)
from RRG_TransCluster import main as program  # noqa: E402
from RRG_TransCluster_Serializer import main as serializer  # noqa: E402

pytestmark = pytest.mark.integration

N, K = 60, 3
N2 = N + 2  # a second size
N_REPS = 6
PROGRAM_ARGS = ["--p-linsp", "(0, 0.5, 6)", "-na", str(N_REPS)]
TILING = ["--ctrl-chunk", "2", "--rep-chunk", "2"]  # 3 x 3 = 9 tiles
N_TILES = 9
TYPED = [str(N), str(K), *PROGRAM_ARGS, *TILING]  # the options of N
SERIAL = ["--N-list", str(N), "--k-list", str(K), *PROGRAM_ARGS, *TILING]
BOTH_SIZES = ["--N-list", str(N), str(N2), *SERIAL[2:]]
TIME_FIELDS = ("t_eig_med", "t_eig_max")
SUBMISSION_KEYS = {
    "utc",
    "job_id",
    "jobname",
    "array",
    "offset",
    "n_tasks",
    "grid",
    "layout",
    "rundir",
    "git",
    "options",
    "command",
}
FIRST_JOB = 4241
# Job ids FIRST_JOB, FIRST_JOB + 1, ...; the n-th call and the later ones
# fail if FAKE_SBATCH_FAIL_FROM=n is set.
FAKE_SBATCH = f"""#!/bin/bash
n=$(( $(cat "$FAKE_SBATCH_COUNT" 2>/dev/null || echo 0) + 1 ))
echo "$n" > "$FAKE_SBATCH_COUNT"
if [ -n "$FAKE_SBATCH_FAIL_FROM" ] && [ "$n" -ge "$FAKE_SBATCH_FAIL_FROM" ]
then
    echo "sbatch: error: boom" >&2
    exit 1
fi
echo "$@" >> "$FAKE_SBATCH_LOG"
cat >> "$FAKE_SBATCH_LOG"
echo "$(( {FIRST_JOB - 1} + n ));clusterA"
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _serializer(root: Path, monkeypatch, *extra: str, argv=SERIAL) -> int:
    monkeypatch.setattr(glue, "PATHDATA", root)
    return serializer([*argv, *extra])


def _print_run(
    root: Path, monkeypatch, capsys, *extra: str
) -> tuple[list[str], str]:
    """The scripts ``--print`` renders, exactly (print adds a newline), and
    its stderr."""
    assert _serializer(root, monkeypatch, "--print", *extra) == TC_EXIT_OK
    out, err = capsys.readouterr()
    bodies = out.split("#!/bin/bash\n")[1:]
    assert all(body.endswith("\n\n") for body in bodies)
    return [f"#!/bin/bash\n{body[:-1]}" for body in bodies], err


def _scripts(root: Path, monkeypatch, capsys, *extra: str) -> list[str]:
    return _print_run(root, monkeypatch, capsys, *extra)[0]


def _status(root: Path, monkeypatch, capsys, *extra: str) -> str:
    assert _serializer(root, monkeypatch, "--status", *extra) == TC_EXIT_OK
    return capsys.readouterr().out


def _expand(spec: str) -> list[int]:
    """Array indices of an ``#SBATCH --array`` value (throttle dropped)."""
    indices = []
    for part in spec.split("%")[0].split(","):
        lo, _, hi = part.partition("-")
        indices += range(int(lo), int(hi or lo) + 1)
    return indices


def _array(script: str) -> tuple[str, str, int]:
    """``(job name, --array value, tile offset)`` of a printed script."""
    name, spec = (
        re.search(rf"^#SBATCH --{key}=(\S+)$", script, re.M)[1]
        for key in ("job-name", "array")
    )
    offset = re.search(r"SLURM_ARRAY_TASK_ID \+ (\d+)\)\)", script)[1]
    return name, spec, int(offset)


def _tasks(script: str) -> list[list[str]]:
    """The program arguments of every task of a printed array script."""
    _, spec, offset = _array(script)
    exec_line = next(ln for ln in script.splitlines() if ln.startswith("exec"))
    words = shlex.split(exec_line)
    assert words[-2] == ser.TILE_FLAG
    return [[*words[3:-2], words[-2], str(offset + i)] for i in _expand(spec)]


def _rundir(root: Path, n_nodes: int = N) -> Path:
    (rundir,) = root.glob(f"*/phtra/N={n_nodes}_k={K}/*")
    return rundir


def _tile_ids(scripts: list[str]) -> list[int]:
    return sorted(int(task[-1]) for s in scripts for task in _tasks(s))


def _tile_paths(tile_id: int) -> TilePaths:
    """The files of tile ``tile_id`` of N (under ``glue.PATHDATA``)."""
    (plan,) = ser.plan_runs([N], [K], TYPED[2:])
    run = plan.run
    tile = decode_tile(run.layout, tile_id)
    return tile_paths(run.rundir, grid_hash(run.grid), tile)


def _assert_agg_equal(a: dict, b: dict) -> None:
    assert sorted(a) == sorted(b)
    for name in ("summary", "meta"):
        for field in a[name].dtype.names:
            if field not in TIME_FIELDS:
                np.testing.assert_array_equal(a[name][field], b[name][field])
    for name in ("ps_ptr", "ps_size", "ps_count", "identity", "cfg"):
        np.testing.assert_array_equal(a[name], b[name])


def _fake_cluster(tmp_path: Path, monkeypatch, fail_from=None) -> Path:
    """``FAKE_SBATCH`` on PATH and a library root in ``tmp_path``."""
    bin_dir, lib = tmp_path / "bin", tmp_path / "lib"
    bin_dir.mkdir()
    lib.mkdir()
    (bin_dir / "sbatch").write_text(FAKE_SBATCH)
    (bin_dir / "sbatch").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_SBATCH_LOG", str(tmp_path / "sbatch.log"))
    monkeypatch.setenv("FAKE_SBATCH_COUNT", str(tmp_path / "sbatch.count"))
    if fail_from is not None:
        monkeypatch.setenv("FAKE_SBATCH_FAIL_FROM", str(fail_from))
    monkeypatch.setattr(ser, "LIB_ROOT", lib)
    monkeypatch.setattr(ser, "PROGRAM", lib / "src" / "RRG_TransCluster.py")
    return lib


# ---------------------------------------------------------------------------
# Arrays, status and resubmission
# ---------------------------------------------------------------------------
def test_array_tasks_reproduce_the_single_process_run(
    tmp_path, monkeypatch, capsys
):
    single, arrays = tmp_path / "single", tmp_path / "arrays"
    monkeypatch.setattr(glue, "PATHDATA", single)
    assert program(TYPED) == TC_EXIT_OK

    scripts = _scripts(arrays, monkeypatch, capsys, "--array-max", "4")
    assert len(scripts) == 3  # tiles 0-3, 4-7 and 8
    tasks = [task for script in scripts for task in _tasks(script)]
    assert sorted(int(t[-1]) for t in tasks) == list(range(N_TILES))
    random.Random(0).shuffle(tasks)
    for task in tasks:
        assert program(task) == TC_EXIT_OK

    status = _status(arrays, monkeypatch, capsys)
    assert f"{N_TILES} complete, 0 partial, 0 missing" in status
    for root in (single, arrays):
        assert _serializer(root, monkeypatch, "--collect", "--boot", "50") == 0
    out, err = capsys.readouterr()
    assert "chi peak" in out
    note = f"no peak interval: {N_REPS} paired reps (need {TC_BOOT_MIN_REPS})"
    assert err.count(note) == 2  # too few reps for a bootstrap, in both runs
    _assert_agg_equal(load_agg(_rundir(single)), load_agg(_rundir(arrays)))


def test_only_missing_resubmits_partial_and_missing_tiles(
    tmp_path, monkeypatch, capsys
):
    everything = _tile_ids(_scripts(tmp_path, monkeypatch, capsys))
    assert everything == list(range(N_TILES))
    for tile_id in (0, 4):
        assert program([*TYPED, "-tid", str(tile_id)]) == TC_EXIT_OK
    _tile_paths(5).part.touch()  # only the name counts

    status = _status(tmp_path, monkeypatch, capsys)
    assert "2 complete, 1 partial, 6 missing" in status
    missing = _tile_ids(
        _scripts(tmp_path, monkeypatch, capsys, "--only-missing")
    )
    assert missing == [1, 2, 3, 5, 6, 7, 8]


def test_only_missing_with_nothing_left_prints_nothing(
    tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    grid = ["--p-list", "0.1", "-na", "2"]
    assert program([str(N), str(K), *grid]) == TC_EXIT_OK
    argv = ["--N-list", str(N), "--k-list", str(K), *grid, "--only-missing"]
    assert serializer(argv) == TC_EXIT_OK
    out, err = capsys.readouterr()
    assert out == "" and "all tiles complete" in err


def test_verify_reads_the_complete_tiles_as_a_task_would(
    tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    assert program(TYPED) == TC_EXIT_OK
    empty, foreign = _tile_paths(0).final, _tile_paths(1).final
    empty.write_bytes(b"")
    write_npz(foreign, {**read_tile(foreign), "identity": np.str_("other")})

    by_name = f"{N_TILES} complete, 0 partial, 0 missing ("
    assert by_name in _status(tmp_path, monkeypatch, capsys)
    verified = "{} complete, 0 partial, 0 missing, {} unreadable, 1 foreign ("
    status = _status(tmp_path, monkeypatch, capsys, "--verify")
    assert verified.format(N_TILES - 2, 1) in status

    scripts, err = _print_run(tmp_path, monkeypatch, capsys, "--only-missing")
    assert scripts == [] and "all tiles complete" in err
    scripts, err = _print_run(
        tmp_path, monkeypatch, capsys, "--only-missing", "--verify"
    )
    assert _tile_ids(scripts) == [0]  # the foreign tile is left alone
    assert f"N={N} k={K}: 1 tiles of another run are not" in err
    (task,) = _tasks(scripts[0])
    with pytest.warns(UserWarning, match="unreadable"):
        assert program(task) == TC_EXIT_OK
    assert len(list(empty.parent.glob(f"{empty.name}.*{TC_BAD_EXT}"))) == 1
    status = _status(tmp_path, monkeypatch, capsys, "--verify")
    assert verified.format(N_TILES - 1, 0) in status


# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------
def test_exec_publishes_submits_and_logs(tmp_path, monkeypatch, capsys):
    lib, data = _fake_cluster(tmp_path, monkeypatch), tmp_path / "d"
    monkeypatch.setenv("SBATCH_TIMELIMIT", "1:00")  # would beat #SBATCH
    scripts = _scripts(data, monkeypatch, capsys, "--array-max", "4")
    arrays = [_array(script) for script in scripts]
    jobs = [str(FIRST_JOB + i) for i in range(len(scripts))]
    assert not data.exists()
    argv = [*SERIAL, "--exec", "--array-max", "4"]
    assert _serializer(data, monkeypatch, argv=argv) == TC_EXIT_OK
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        f"{name}: job {job}, {len(_expand(spec))} tasks"
        for job, (name, spec, _) in zip(jobs, arrays)
    ]
    assert "SBATCH_TIMELIMIT" in err and "override the #SBATCH lines" in err
    assert (lib / SBATCH_LOG_DIR).is_dir()
    submitted = (tmp_path / "sbatch.log").read_text()
    assert submitted == "".join(f"--parsable\n{script}" for script in scripts)

    (plan,) = ser.plan_runs([N], [K], TYPED[2:])
    rundir = _rundir(data)
    assert rundir == plan.run.rundir and (rundir / CFG_SIDECAR_NAME).is_file()
    logs = sorted((rundir / TC_SUBMISSIONS_DIR).iterdir())
    assert [path.name for path in logs] == [f"{job}.json" for job in jobs]
    for path, job, (name, spec, offset) in zip(logs, jobs, arrays):
        record = json.loads(path.read_text())
        assert set(record) == SUBMISSION_KEYS
        assert (record["job_id"], record["jobname"]) == (job, name)
        assert (record["array"], record["offset"]) == (spec, offset)
        assert record["n_tasks"] == len(_expand(spec))
        assert record["layout"] == dataclasses.asdict(plan.run.layout)
        assert record["grid"] == grid_hash(plan.run.grid)
        assert record["rundir"] == str(rundir)
        git = record["git"]  # None outside a git checkout
        assert git is None or sorted(git) == ["commit", "dirty"]
        assert record["options"][: len(TYPED)] == TYPED
        assert record["options"] == list(plan.options)
        assert record["command"][1:] == argv


@pytest.mark.parametrize("fail_from", [1, 2])
def test_exec_lists_the_arrays_a_failed_sbatch_left_out(
    fail_from, tmp_path, monkeypatch, capsys
):
    _fake_cluster(tmp_path, monkeypatch, fail_from)
    data = tmp_path / "d"
    scripts = _scripts(data, monkeypatch, capsys, "--array-max", "4")
    code = _serializer(data, monkeypatch, "--exec", "--array-max", "4")
    assert code == TC_EXIT_FAILED
    out, err = capsys.readouterr()
    left_out = [
        f"  N={N} k={K} {name}: tiles {offset} + [{spec}]"
        for name, spec, offset in map(_array, scripts[fail_from - 1 :])
    ]
    assert "boom" in err and err.splitlines()[-len(left_out) :] == left_out
    assert (tmp_path / "sbatch.count").read_text().strip() == str(fail_from)
    logged = (_rundir(data) / TC_SUBMISSIONS_DIR).glob("*.json")
    jobs = [str(FIRST_JOB + i) for i in range(fail_from - 1)]
    assert sorted(path.stem for path in logged) == jobs
    assert len(out.splitlines()) == len(jobs)


@pytest.mark.parametrize("fail_at", [2, 3])  # a middle and the last array
def test_exec_stops_at_a_submission_log_it_cannot_write(
    fail_at, tmp_path, monkeypatch, capsys
):
    _fake_cluster(tmp_path, monkeypatch)
    data = tmp_path / "d"
    scripts = _scripts(data, monkeypatch, capsys, "--array-max", "4")
    arrays = [_array(script) for script in scripts]
    log, jobs = ser.log_submission, []

    def log_until_the_quota_is_full(plan, array, job_id, command):
        jobs.append(job_id)
        if len(jobs) == fail_at:
            raise OSError(errno.EDQUOT, os.strerror(errno.EDQUOT))
        return log(plan, array, job_id, command)

    monkeypatch.setattr(ser, "log_submission", log_until_the_quota_is_full)
    code = _serializer(data, monkeypatch, "--exec", "--array-max", "4")
    assert code == TC_EXIT_FAILED
    out, err = capsys.readouterr()
    assert out.splitlines() == [  # the job whose log failed too
        f"{name}: job {job}, {len(_expand(spec))} tasks"
        for job, (name, spec, _) in zip(jobs, arrays)
    ]
    assert f"error: job {jobs[-1]} is queued, its log was not written" in err
    left_out = [
        f"  N={N} k={K} {name}: tiles {offset} + [{spec}]"
        for name, spec, offset in arrays[fail_at:]
    ]
    listed = [line for line in err.splitlines() if line.startswith("  N=")]
    assert listed == left_out and ("not submitted" in err) == bool(left_out)
    assert (tmp_path / "sbatch.count").read_text().strip() == str(fail_at)
    logged = (_rundir(data) / TC_SUBMISSIONS_DIR).glob("*.json")
    assert sorted(path.stem for path in logged) == jobs[:-1]


def test_exec_submits_nothing_if_a_run_directory_holds_another_run(
    tmp_path, monkeypatch, capsys
):
    _fake_cluster(tmp_path, monkeypatch)
    data = tmp_path / "d"
    monkeypatch.setattr(glue, "PATHDATA", data)
    _, second = ser.plan_runs([N, N2], [K], TYPED[2:])
    publish_cfg_sidecar(second.run.rundir, graph={}, model={}, run={"x": 1})
    code = _serializer(data, monkeypatch, "--exec", argv=BOTH_SIZES)
    assert code == TC_EXIT_CFG_MISMATCH
    assert "exists with different content" in capsys.readouterr().err
    assert not (tmp_path / "sbatch.count").exists()  # sbatch never ran
    assert not any(data.rglob(TC_SUBMISSIONS_DIR))


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores permission bits")
def test_exec_submits_nothing_if_a_run_directory_cannot_take_a_log(
    tmp_path, monkeypatch, capsys
):
    _fake_cluster(tmp_path, monkeypatch)
    data = tmp_path / "d"
    monkeypatch.setattr(glue, "PATHDATA", data)
    _, second = ser.plan_runs([N, N2], [K], TYPED[2:])
    logs = second.run.rundir / TC_SUBMISSIONS_DIR
    logs.mkdir(parents=True)
    logs.chmod(0o500)  # read-only
    try:
        code = _serializer(data, monkeypatch, "--exec", argv=BOTH_SIZES)
    finally:
        logs.chmod(0o700)
    assert code == TC_EXIT_FAILED
    err = capsys.readouterr().err
    assert f"error: [Errno {errno.EACCES}] " in err and str(logs) in err
    assert not (tmp_path / "sbatch.count").exists()  # not even for N


def test_exec_without_sbatch_submits_nothing(tmp_path, monkeypatch, capsys):
    lib, data, empty = tmp_path / "lib", tmp_path / "d", tmp_path / "bin"
    lib.mkdir()
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))  # never a real sbatch
    monkeypatch.setattr(ser, "LIB_ROOT", lib)
    monkeypatch.setattr(ser, "PROGRAM", lib / "src" / "RRG_TransCluster.py")
    assert _serializer(data, monkeypatch, "--exec") == TC_EXIT_FAILED
    assert "error: sbatch not found" in capsys.readouterr().err
    assert not data.exists()  # no cfg.json
    assert not (lib / SBATCH_LOG_DIR).exists()


# ---------------------------------------------------------------------------
# Collect
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "damage, message", [("empty", "unreadable"), ("clash", "disagree on n_neg")]
)
def test_collect_names_a_bad_tile_and_goes_on(
    damage, message, tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    for n_nodes in (N, N2):
        assert program([str(n_nodes), *TYPED[1:]]) == TC_EXIT_OK
    tile = _tile_paths(0).final
    if damage == "empty":
        tile.write_bytes(b"")
        named = [tile.name]
    else:  # a copy that disagrees on the units it shares with the tile
        data = read_tile(tile)
        records = data["records"].copy()
        records["n_neg"] += 1
        copy = tile.with_name(f"copy_{tile.name}")
        write_npz(copy, {**data, "records": records})
        named = [tile.name, copy.name]

    code = _serializer(tmp_path, monkeypatch, "--collect", argv=BOTH_SIZES)
    assert code == TC_EXIT_FAILED
    out, err = capsys.readouterr()
    (line,) = [ln for ln in err.splitlines() if ln.startswith(f"N={N} k=")]
    assert message in line and all(name in line for name in named)
    assert not (_rundir(tmp_path) / TC_AGG_NAME).exists()
    agg = _rundir(tmp_path, N2) / TC_AGG_NAME
    assert agg.is_file() and f"N={N2} k={K}: {agg}" in out.splitlines()


def test_collect_reports_runs_without_tiles(tmp_path, monkeypatch, capsys):
    code = _serializer(tmp_path, monkeypatch, "--collect")
    assert code == TC_EXIT_FAILED
    assert f"N={N} k={K}" in capsys.readouterr().err


@pytest.mark.parametrize(
    "tile_ids, window, reason",
    [
        (  # p = 0.2 and 0.3 only
            range(N_TILES),
            ["--peak-window", "0.15", "0.35"],
            "2 control values with converged reps in the window (need 3)",
        ),
        (  # reps 0-1 at p <= 0.3, reps 2-3 at p >= 0.4: none has them all
            (0, 1, 5),
            [],
            "chi is flat or undefined over the window (0 paired reps)",
        ),
    ],
)
def test_collect_says_why_there_is_no_chi_peak(
    tile_ids, window, reason, tmp_path, monkeypatch, capsys
):
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    for tile_id in tile_ids:
        assert program([*TYPED, "-tid", str(tile_id)]) == TC_EXIT_OK
    code = _serializer(tmp_path, monkeypatch, "--collect", *window)
    assert code == TC_EXIT_OK  # agg.npz is valid all the same
    out, err = capsys.readouterr()
    assert (_rundir(tmp_path) / TC_AGG_NAME).is_file()
    assert "chi peak nan [nan, nan]" in out
    assert f"N={N} k={K}: no chi peak: {reason}" in err.splitlines()
