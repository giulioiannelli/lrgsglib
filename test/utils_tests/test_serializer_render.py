"""SLURM array rendering (kernels/Serializer) and RRG serializer planning."""

import math
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[2] / "src"
sys.path.insert(0, str(SRC))

import kernels.RRG_TransCluster as glue  # noqa: E402
import kernels.RRG_TransCluster_Serializer as ser  # noqa: E402
from kernels.Serializer import (  # noqa: E402
    check_sbatch_value,
    format_array_spec,
    render_sbatch_tiles_script,
    split_arrays,
)
from kernels.TransClusterStats import META_DTYPE, SUMMARY_DTYPE  # noqa: E402
from lrgsglib.config.progargs.defs.generic import (  # noqa: E402
    SBATCH_FORBIDDEN_CHARS,
)
from lrgsglib.config.progargs.defs.RandomRegular import (  # noqa: E402
    DEFAULT_RRG_NAVG,
    DEFAULT_RRG_SEED,
    DEFAULT_RRG_WORKDIR,
    RRG_SAMPLER_NATIVE,
    RRG_SAMPLER_PAIRING,
)
from lrgsglib.config.progargs.defs.TransCluster import (  # noqa: E402
    DEFAULT_TC_BOOT,
    DEFAULT_TC_CTRL_CHUNK,
    DEFAULT_TC_GAP_EVERY,
    DEFAULT_TC_REP_CHUNK,
    DEFAULT_TC_REP_START,
    TC_BOOT_MIN_REPS,
    TC_BOOT_QUANTILES,
    TC_BOOT_SEED_LIMIT,
    TC_MEM_BASE_MB,
    TC_MEM_ENGINE_MB_PER_NODE,
    TC_MEM_MB_PER_NODE,
    TC_MEM_SAFETY,
    TC_PEAK_MIN_CTRLS,
    TC_SPLIT_CTRL_CHUNK,
    TC_SPLIT_CTRL_MIN_N,
    TC_SRUN_ACTIONS,
)
from RRG_TransCluster_Serializer import main  # noqa: E402

BASH = shutil.which("bash")
WORDS = ["a b", "semi;colon", "gt>lt<", "it's", "$HOME", "`id`", "-x", ""]
RENDER = {
    "jobname": "RRGTC_N60_k3_a0",
    "chdir": "/abs/lib",
    "array": "0-3,7%2",
    "mem_mb": 700,
    "cpus": 1,
    "time": "01:00:00",
    "partition": "",
    "signal": "B:USR1@300",
    "env": {"OMP_NUM_THREADS": "1", "LRGSG_DATA": "/data dir/it's"},
    "command": ["python", "src/P.py", "60", "3"],
    "index_flag": "--tile-id",
    "tile_offset": 7,
    "log_dir": ".log",
}
GOLDEN = """#!/bin/bash
#SBATCH --job-name=RRGTC_N60_k3_a0
#SBATCH --chdir=/abs/lib
#SBATCH --output=.log/%x_%A_%a.out
#SBATCH --error=.log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --array=0-3,7%2
#SBATCH --mem=700M
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --signal=B:USR1@300

set -euo pipefail
export OMP_NUM_THREADS=1 LRGSG_DATA='/data dir/it'"'"'s'
exec python src/P.py 60 3 --tile-id "$((SLURM_ARRAY_TASK_ID + 7))"
"""
ENGINE = ["-ge", "nx"]  # planning then needs no graph_tool
BASE = ["--N-list", "60", "--k-list", "3", "--p-linsp", "(0, 0.5, 6)", *ENGINE]
GRID = ["--p-list", "0.1", "0.2", *ENGINE]
# Values an #SBATCH line would not pass on verbatim: split at whitespace,
# ended at #, quotes and backslashes stripped; the required ones blank.
BAD_SBATCH = [
    ("jobname", "a b"),
    ("chdir", "/x y"),
    ("log_dir", ""),
    ("array", ""),
    ("signal", "B:USR1@300 --qos=x"),
    ("time", "12:00 --qos=x"),
    ("partition", "p\tq"),
    *(("jobname", f"a{char}b") for char in SBATCH_FORBIDDEN_CHARS),
    *(("partition", f"p{char}") for char in SBATCH_FORBIDDEN_CHARS),
]
# A new code default for each PINNED option, as if changed between exec and
# the start of a task. Identity options (--eig-tol, -ge) are left out: a
# changed default of theirs fails on cfg.json instead.
LATER_DEFAULTS = {
    "number_of_averages": DEFAULT_RRG_NAVG + 1,
    "rep_start": DEFAULT_TC_REP_START + 1,
    "ctrl_chunk": DEFAULT_TC_CTRL_CHUNK + 1,
    "rep_chunk": DEFAULT_TC_REP_CHUNK + 1,
    "seed": DEFAULT_RRG_SEED + 1,
    "workdir": f"{DEFAULT_RRG_WORKDIR}_later",
    "sampler": RRG_SAMPLER_NATIVE,
    "gap_every": DEFAULT_TC_GAP_EVERY + 1,
}
# The fewest resamples at a window end that warn: the tail mass beyond either
# TC_BOOT_QUANTILES, from which the chi-peak interval may be clipped.
N_EDGE = math.ceil(
    min(TC_BOOT_QUANTILES[0], 1 - TC_BOOT_QUANTILES[1]) * DEFAULT_TC_BOOT
)
BOUND = "the chi maximum is at a window end; the location is only a bound"
CLIPPED = (
    f"{N_EDGE / DEFAULT_TC_BOOT:.0%} of the usable resamples peak at a "
    "window end; the interval may be clipped"
)
NO_BOOT = {  # the META fields of --boot 0
    "peak_lo": np.nan,
    "peak_hi": np.nan,
    "peak_sd": np.nan,
    "peak_edge_frac": np.nan,
    "n_boot": 0,
}


# ---------------------------------------------------------------------------
# Array specs and splitting
# ---------------------------------------------------------------------------
def test_array_spec_compresses_runs_and_caps_the_throttle():
    ids = [*range(10), 17, 42, 43]
    assert format_array_spec(ids, 5) == "0-9,17,42-43%5"
    assert format_array_spec([43, 42, 42], 99) == "42-43%2"
    assert format_array_spec([7], 3) == "7%1"


@pytest.mark.parametrize(
    "ids, concurrent", [([], 1), ([-1, 0], 1), ([0, 1], 0)]
)
def test_array_spec_rejects_bad_input(ids, concurrent):
    with pytest.raises(ValueError):
        format_array_spec(ids, concurrent)


@pytest.mark.parametrize(
    "ids, array_max",
    [(range(2500), 1000), ([5, 999, 1000, 2999, 5], 1000), (range(7), 3)],
)
def test_split_arrays_covers_every_tile_below_the_index_cap(ids, array_max):
    arrays = split_arrays(ids, array_max)
    covered = [off + i for off, indices in arrays for i in indices]
    assert covered == sorted(set(ids))
    for offset, indices in arrays:
        assert offset % array_max == 0 and indices
        assert all(0 <= i < array_max for i in indices)


def test_split_arrays_offsets():
    assert split_arrays([5, 999, 1000, 2999], 1000) == [
        (0, [5, 999]),
        (1000, [0]),
        (2000, [999]),
    ]
    with pytest.raises(ValueError):
        split_arrays([0], 0)
    with pytest.raises(ValueError):
        split_arrays([-1], 10)


# ---------------------------------------------------------------------------
# Script rendering
# ---------------------------------------------------------------------------
def test_render_golden_script():
    assert render_sbatch_tiles_script(**RENDER) == GOLDEN


def test_render_optional_lines():
    script = render_sbatch_tiles_script(
        **{**RENDER, "time": "", "partition": "p"}
    )
    assert "--time" not in script and "#SBATCH --partition=p\n" in script


@pytest.mark.parametrize("field, value", BAD_SBATCH)
def test_render_rejects_values_sbatch_would_not_pass_verbatim(field, value):
    with pytest.raises(ValueError, match=field):
        render_sbatch_tiles_script(**{**RENDER, field: value})


def test_optional_sbatch_values_may_be_blank():
    check_sbatch_value("time", "", required=False)
    with pytest.raises(ValueError, match="time must be non-empty"):
        check_sbatch_value("time", "")


def _run_script(tmp_path: Path, env_value: str, task: str | None):
    """Render a script whose task prints its arguments, then run it."""
    printer = 'printf "%s\\n" "$LRGSG_DATA" "$@"'
    script = render_sbatch_tiles_script(
        **{
            **RENDER,
            "env": {"LRGSG_DATA": env_value},
            "command": ["sh", "-c", printer, "sh", *WORDS],
        }
    )
    path = tmp_path / "task.sh"
    path.write_text(script)
    env = {"PATH": "/usr/bin:/bin"}
    if task is not None:
        env["SLURM_ARRAY_TASK_ID"] = task
    return subprocess.run(
        [BASH, str(path)], env=env, capture_output=True, text=True
    )


@pytest.mark.skipif(BASH is None, reason="needs bash")
def test_rendered_script_passes_every_word_through_the_shell(tmp_path):
    value = "/data dir/it's $x;"
    result = _run_script(tmp_path, value, task="3")
    assert result.returncode == 0, result.stderr
    assert result.stdout.split("\n")[:-1] == [value, *WORDS, "--tile-id", "10"]


@pytest.mark.skipif(BASH is None, reason="needs bash")
def test_rendered_script_fails_outside_an_array(tmp_path):
    result = _run_script(tmp_path, "/data", task=None)
    assert result.returncode != 0 and result.stdout == ""


# ---------------------------------------------------------------------------
# RRG serializer planning
# ---------------------------------------------------------------------------
def test_actions_match_the_command_line_vocabulary():
    assert set(ser.ACTIONS) == set(TC_SRUN_ACTIONS)


def _task_args(plan: ser.Plan):
    """The arguments every task of ``plan`` runs with (tile id aside)."""
    return ser.program_parser().parse_args(list(plan.options))


@pytest.mark.parametrize(
    "extra, ctrl_chunk",
    [
        ([], TC_SPLIT_CTRL_CHUNK),
        (["--ctrl-chunk", "0"], 0),
        (["--ctrl-chunk=3"], 3),
    ],
)
def test_large_runs_split_over_the_control_value(extra, ctrl_chunk):
    plan = ser.plan_run(TC_SPLIT_CTRL_MIN_N, 3, [*GRID, *extra])
    assert _task_args(plan).ctrl_chunk == ctrl_chunk
    assert sum(w.startswith("--ctrl-chunk") for w in plan.options) == 1
    small = ser.plan_run(TC_SPLIT_CTRL_MIN_N - 2, 3, GRID)
    assert _task_args(small).ctrl_chunk == DEFAULT_TC_CTRL_CHUNK


def test_tasks_pin_the_resolved_values_of_untyped_options():
    """A task that starts later keeps the planned values even if a code
    default changes meanwhile. Typed options stay as typed, once."""
    typed = [*GRID, "--seed", "5", "--rep-chunk=2"]
    plan = ser.plan_run(60, 3, typed)
    head, tail = plan.options[: 2 + len(typed)], plan.options[2 + len(typed) :]
    assert head == ("60", "3", *typed)
    now = vars(_task_args(plan))
    untyped = [d for d in ser.PINNED if d not in ("seed", "rep_chunk")]
    assert tail[::2] == tuple(ser.PROGRAM_FLAGS[d] for d in untyped)
    assert tail[1::2] == tuple(str(now[d]) for d in untyped)
    assert now["sampler"] == plan.run.method == RRG_SAMPLER_PAIRING
    later = ser.program_parser()
    later.set_defaults(**{dest: object() for dest in ser.PINNED})
    changed = vars(later.parse_args(list(plan.options)))
    assert [d for d in ser.PINNED if changed[d] != now[d]] == []
    assert plan.options.count("--seed") == 1
    # identity options stay as typed: a changed default fails on cfg.json
    identity = {"eig_tol", "graph_engine"}  # --eig-tol, -ge
    assert identity <= set(ser.PROGRAM_FLAGS) and not identity & set(ser.PINNED)


def test_a_later_code_default_leaves_queued_tasks_on_their_tiles():
    """A task that starts after a code default changed keeps the units of
    its tile id, its seed, run directory, sampler and lambda_1 subsample."""
    assert set(ser.PINNED) >= set(LATER_DEFAULTS)
    plan = ser.plan_run(60, 3, GRID)  # none of them typed
    later = ser.program_parser()
    later.set_defaults(**LATER_DEFAULTS)
    argv = [*plan.options, ser.TILE_FLAG, "0"]
    run = glue.resolve_run(later.parse_args(argv), [str(ser.PROGRAM), *argv])
    for field in ("layout", "seed", "rundir", "method", "gap_every"):
        assert getattr(run, field) == getattr(plan.run, field), field


def test_forwarded_options_must_be_program_options():
    ser.check_forwarded(["--p-list", "-1e-3", "--seed=3", "(0, 1, 2)", "-v"])
    with pytest.raises(ValueError, match=r"^unrecognized option --qqqqqqqq$"):
        ser.check_forwarded(["--qqqqqqqq"])


def test_boolean_program_options_are_forwarded():
    # -v is a BooleanOptionalAction: the program also accepts --no-verbose.
    ser.check_forwarded(["--no-verbose", "--verbose"])


@pytest.mark.parametrize(
    "flag", [["-tid", "3"], ["--tile-id", "3"], ["--tile-id=0"], ["-tid", "-1"]]
)
def test_the_tile_id_belongs_to_the_serializer(flag, capsys):
    with pytest.raises(SystemExit) as exit_info:
        main([*BASE, *flag])
    assert exit_info.value.code == 2
    assert "set per task" in capsys.readouterr().err


@pytest.mark.parametrize(
    "extra, message",
    [
        (["--N-list", "0"], "--N-list and --k-list take positive integers"),
        (["--array-max", "0"], "--array-max and --array-concurrent must be"),
        (["-ac", "0"], "--array-max and --array-concurrent must be >= 1"),
        (["--mem-mb", "-1"], "--mem-mb and --boot must be >= 0"),
        (["--collect", "--boot", "-1"], "--mem-mb and --boot must be >= 0"),
        (
            ["--collect", "--boot-seed", "-1"],
            f"--boot-seed must lie in [0, {TC_BOOT_SEED_LIMIT})",
        ),
        (
            ["--collect", "--boot-seed", str(TC_BOOT_SEED_LIMIT)],
            "--boot-seed must lie in",
        ),
        (["--collect", "--peak-window", "0.3", "0.1"], "needs LO < HI"),
        (["--collect", "--peak-window", "0.3", "0.3"], "needs LO < HI"),
        (["--collect", "--peak-window", "nan", "0.3"], "needs LO < HI"),
        (["--collect", "--peak-window", "0.3", "nan"], "needs LO < HI"),
        (["--status", "--boot", "5"], "--boot does not apply to --status"),
        (["--status", "-ac", "5"], "--array-concurrent does not apply to"),
        (["--collect", "--only-missing"], "--only-missing does not apply"),
        (["--exec", "--legacy-txt"], "--legacy-txt does not apply to --exec"),
        (["--verify"], "--verify applies to --status and --only-missing"),
        (["--collect", "--verify"], "--verify applies to --status and"),
        (
            ["-aT", "12:00 --qos=x"],
            "--array-time '12:00 --qos=x' goes into an #SBATCH line",
        ),
        (["-ap", "a#b"], "--array-partition 'a#b' goes into an #SBATCH line"),
        (
            ["--arry-max", "4"],
            "unrecognized option --arry-max (did you mean --array-max?)",
        ),
        (["--sed=5"], "unrecognized option --sed (did you mean --seed?)"),
    ],
)
def test_serializer_options_are_checked(extra, message, capsys):
    with pytest.raises(SystemExit) as exit_info:
        main([*BASE, *extra])
    assert exit_info.value.code == 2
    assert message in capsys.readouterr().err


def test_a_library_path_sbatch_would_split_is_refused(monkeypatch, capsys):
    monkeypatch.setattr(ser, "LIB_ROOT", Path("/lib dir"))
    with pytest.raises(SystemExit) as exit_info:
        main([*BASE, "--print"])
    assert exit_info.value.code == 2
    assert "the library path '/lib dir' goes into" in capsys.readouterr().err


def test_help_lists_the_default_action_once_and_no_positionals(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["-h"])
    assert exit_info.value.code == 0
    out = capsys.readouterr().out
    assert "(default: print)" not in out and "positional arguments" not in out
    assert all(f"--{action}" in out for action in TC_SRUN_ACTIONS)
    assert ser.TILE_FLAG not in out  # set per task


def test_task_memory_estimate():
    per_node = TC_MEM_MB_PER_NODE + TC_MEM_ENGINE_MB_PER_NODE["gt"]
    expected = math.ceil(TC_MEM_SAFETY * (TC_MEM_BASE_MB + 1000 * per_node))
    assert ser.task_mem_mb(1000, "gt", 0) == expected
    assert ser.task_mem_mb(10**6, "nx", 0) > ser.task_mem_mb(10**6, "gt", 0)
    assert ser.task_mem_mb(1000, "gt", 1234) == 1234


def test_print_is_side_effect_free(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(glue, "PATHDATA", tmp_path / "data")
    log_dir = ser.LIB_ROOT / ".log"
    had_log_dir = log_dir.exists()
    assert main([*BASE, "-na", "4", "--rep-chunk", "1", "--print"]) == 0
    out = capsys.readouterr().out
    assert not (tmp_path / "data").exists()
    assert log_dir.exists() == had_log_dir
    (script,) = out.split("#!/bin/bash\n")[1:]
    assert "#SBATCH --array=0-3%4\n" in script
    exec_line = next(ln for ln in script.splitlines() if ln.startswith("exec"))
    words = shlex.split(exec_line)
    assert words[1:3] == [sys.executable, "src/RRG_TransCluster.py"]
    assert words[3:] == [
        "60",
        "3",
        "--p-linsp",
        "(0, 0.5, 6)",
        *ENGINE,
        "-na",
        "4",
        "--rep-chunk",
        "1",
        # the resolved values of the untyped PINNED options
        "--rep-start",
        str(DEFAULT_TC_REP_START),
        "--ctrl-chunk",
        str(DEFAULT_TC_CTRL_CHUNK),
        "--seed",
        str(DEFAULT_RRG_SEED),
        "--workdir",
        DEFAULT_RRG_WORKDIR,
        "--sampler",
        RRG_SAMPLER_PAIRING,
        "--gap-every",
        str(DEFAULT_TC_GAP_EVERY),
        "--tile-id",
        "$((SLURM_ARRAY_TASK_ID + 0))",
    ]
    data = shlex.quote(str(tmp_path / "data"))
    assert f"LRGSG_DATA={data}" in script


def test_print_names_one_array_per_block(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(glue, "PATHDATA", tmp_path)
    argv = [*BASE, "-na", "5", "--rep-chunk", "1", "--array-max", "2"]
    assert main([*argv, "-ac", "1"]) == 0
    scripts = capsys.readouterr().out.split("#!/bin/bash\n")[1:]
    assert len(scripts) == 3
    for block, script in enumerate(scripts):
        assert f"--job-name=RRGTC_N60_k3_a{block}\n" in script
        assert f"+ {2 * block}))" in script
    assert "#SBATCH --array=0%1\n" in scripts[-1]


# ---------------------------------------------------------------------------
# Chi-peak report
# ---------------------------------------------------------------------------
def _meta(**fields) -> np.ndarray:
    """A META row of ``collect``: an interior chi peak from 30 paired reps
    and ``DEFAULT_TC_BOOT`` resamples, none at a window end, unless
    ``fields`` say otherwise."""
    row = dict.fromkeys(META_DTYPE.names, 0) | {
        "peak": 0.3,
        "peak_lo": 0.29,
        "peak_hi": 0.31,
        "peak_reps": 30,
        "window_lo": -np.inf,
        "window_hi": np.inf,
        "n_boot": DEFAULT_TC_BOOT,
        **fields,
    }
    return np.array(tuple(row[name] for name in META_DTYPE.names), META_DTYPE)


@pytest.mark.parametrize(
    "fields, warning",
    [
        ({}, ""),
        ({"peak_edge_frac": 1 / DEFAULT_TC_BOOT}, ""),  # one resample
        ({"peak_edge_frac": (N_EDGE - 1) / DEFAULT_TC_BOOT}, ""),  # one fewer
        ({"peak_edge_frac": N_EDGE / DEFAULT_TC_BOOT}, CLIPPED),
        ({"peak_at_edge": 1, "peak_edge_frac": 1.0}, BOUND),
        (NO_BOOT, ""),  # no resample, so no resample clause
        ({**NO_BOOT, "peak_at_edge": 1}, BOUND),
    ],
)
def test_the_peak_warns_only_if_a_window_end_bounds_it(fields, warning, capsys):
    ser._print_peak(_meta(**fields))
    out, err = capsys.readouterr()
    level = TC_BOOT_QUANTILES[1] - TC_BOOT_QUANTILES[0]
    assert out.startswith("chi peak 0.3 [")
    assert out.endswith(f"] ({level:.0%}) from 30 paired reps\n")
    hint = " (widen --peak-window or the grid)\n"
    assert err == (f"warning: {warning}{hint}" if warning else "")


@pytest.mark.parametrize(
    "peak_reps, reason",
    [
        (
            TC_BOOT_MIN_REPS - 1,
            f"{TC_BOOT_MIN_REPS - 1} paired reps (need {TC_BOOT_MIN_REPS})",
        ),
        (TC_BOOT_MIN_REPS, "chi is flat in every resample"),
    ],
)
def test_no_interval_reason(peak_reps, reason):
    assert ser._no_interval_reason(_meta(peak_reps=peak_reps)) == reason


def test_no_peak_reason_counts_only_control_values_with_reps():
    """A control value in the window whose reps all failed does not count."""
    summary = np.zeros(TC_PEAK_MIN_CTRLS, SUMMARY_DTYPE)
    summary["ctrl"] = np.linspace(0.1, 0.3, TC_PEAK_MIN_CTRLS)
    summary["n_rep"] = 1
    summary["n_rep"][1] = 0
    assert ser._no_peak_reason(summary, _meta()) == (
        f"{TC_PEAK_MIN_CTRLS - 1} control values with converged reps in the "
        f"window (need {TC_PEAK_MIN_CTRLS})"
    )
