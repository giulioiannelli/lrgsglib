"""RRG_TransCluster_Serializer: SLURM arrays of RRG_TransCluster tiles.

One run per (N, k). Its tiles are split into arrays of at most
``--array-max`` tasks; the task with array index ``i`` runs

    RRG_TransCluster.py N k <options> --tile-id <offset + i>

Options the serializer does not know are forwarded to every task. They are
first checked by the program's own parser and ``resolve_run``, so the tile
layout, run directory and identity used here are the program's. The tasks
also get the resolved value of every ``PINNED`` option that was not typed,
so a task that starts days later cannot pick up a changed code default.
From ``TC_SPLIT_CTRL_MIN_N`` nodes on, tiles hold ``TC_SPLIT_CTRL_CHUNK``
control values unless ``--ctrl-chunk`` is given.

Actions: ``print`` (default, no side effects), ``exec`` (check every
cfg.json and ``submissions/`` directory, then ``sbatch`` and one
``submissions/<job id>.json`` per array),
``status`` (tile counts) and ``collect`` (``agg.npz`` and a cost table).
``--only-missing`` restricts print and exec to the tiles a task would
compute; ``--verify`` reads every complete tile instead of trusting its
name.
"""

from __future__ import annotations

import argparse
import dataclasses
import difflib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import kernels.RRG_TransCluster as glue
from kernels.Serializer import (
    build_jobname,
    check_sbatch_value,
    format_array_spec,
    render_sbatch_tiles_script,
    split_arrays,
)
from kernels.TransClusterStats import collect, cost_table
from kernels.TransClusterTiles import (
    Tile,
    TileError,
    TileHeader,
    TileMismatch,
    TileUnreadable,
    check_tile,
    decode_tile,
    git_state,
    grid_hash,
    read_tile,
    tile_paths,
)
from lrgsglib.config.progargs.defs.generic import (
    DEFAULT_SRUN_THREADS,
    SBATCH_ENV,
    SBATCH_ENV_PREFIX,
    SBATCH_LOG_DIR,
    SRUN_THREAD_ENV_VARS,
)
from lrgsglib.config.progargs.defs.TransCluster import (
    TC_AGG_NAME,
    TC_BOOT_MIN_REPS,
    TC_BOOT_QUANTILES,
    TC_BOOT_SEED_LIMIT,
    TC_EXIT_CFG_MISMATCH,
    TC_EXIT_FAILED,
    TC_EXIT_OK,
    TC_MEM_BASE_MB,
    TC_MEM_ENGINE_MB_PER_NODE,
    TC_MEM_MB_PER_NODE,
    TC_MEM_SAFETY,
    TC_PEAK_MIN_CTRLS,
    TC_SIGNAL_LEAD_SECONDS,
    TC_SPLIT_CTRL_CHUNK,
    TC_SPLIT_CTRL_MIN_N,
    TC_STOP_SIGNALS,
    TC_SUBMISSIONS_DIR,
    TC_TILES_DIR,
)
from lrgsglib.config.progargs.RandomRegular import (
    RRG_TC_TILE_FLAGS,
    RRG_TransCluster_opt_args,
    RRG_TransCluster_progName,
    RRG_TransCluster_progNameShrt,
    RRG_TransCluster_srun_actions,
    RRG_TransCluster_srun_args,
)
from lrgsglib.statsys._io import atomic_write_bytes
from lrgsglib.statsys._naming import publish_cfg_sidecar
from parsers.RRG_TransCluster import build_parser as build_program_parser

__all__ = [
    "ACTIONS",
    "Array",
    "EDGE_FRAC_WARN",
    "LIB_ROOT",
    "OPTION_ACTIONS",
    "PINNED",
    "PROGRAM",
    "PROGRAM_FLAGS",
    "Plan",
    "TILE_FLAG",
    "check_forwarded",
    "check_options",
    "collect_runs",
    "log_submission",
    "pending_tiles",
    "plan_run",
    "plan_runs",
    "print_scripts",
    "program_parser",
    "report_status",
    "sbatch",
    "submit_arrays",
    "task_arrays",
    "task_env",
    "task_mem_mb",
    "task_options",
    "tile_states",
]

SRC_DIR = Path(__file__).resolve().parents[1]
LIB_ROOT = SRC_DIR.parent  # the tasks' --chdir
PROGRAM = SRC_DIR / f"{RRG_TransCluster_progName}.py"
# dest -> long flag of each program option (the long flag comes last)
PROGRAM_FLAGS = {
    flags[-1].lstrip("-").replace("-", "_"): flags[-1]
    for flags in RRG_TransCluster_opt_args
}
TILE_FLAG = RRG_TC_TILE_FLAGS[-1]
# Resolved values every task gets unless typed: the tile id -> unit map
# (reps and chunks), seed, workdir, sampler and lambda_1 subsample. Identity
# values whose defaults are elided from the run directory name (--eig-tol,
# -ge, ...) stay as typed: a changed default then fails loudly on cfg.json
# instead of starting a new run directory.
PINNED = (
    "number_of_averages",
    "rep_start",
    "ctrl_chunk",
    "rep_chunk",
    "seed",
    "workdir",
    "sampler",
    "gap_every",
)
_STOP = TC_STOP_SIGNALS[0].removeprefix("SIG")
STOP_SIGNAL = f"B:{_STOP}@{TC_SIGNAL_LEAD_SECONDS}"  # to the exec'd program
TILE_STATES = ("complete", "partial", "missing")  # judged by file name
VERIFIED_STATES = (*TILE_STATES, "unreadable", "foreign")  # --verify
RECOMPUTED = ("partial", "missing", "unreadable")  # what a task computes
SCRIPT_ACTIONS = ("print", "exec")  # the actions that render sbatch scripts
# The actions that read each serializer option. The others refuse it unless
# it is left at its default, since they would silently ignore it.
OPTION_ACTIONS = {
    **dict.fromkeys(
        (
            "only_missing",
            "array_concurrent",
            "array_partition",
            "array_time",
            "array_max",
            "mem_mb",
        ),
        SCRIPT_ACTIONS,
    ),
    **dict.fromkeys(
        ("boot", "boot_seed", "peak_window", "legacy_txt"), ("collect",)
    ),
}
# Fraction of the resamples peaking at a window end from which the chi-peak
# interval may be clipped: the tail mass beyond either TC_BOOT_QUANTILES.
EDGE_FRAC_WARN = min(TC_BOOT_QUANTILES[0], 1 - TC_BOOT_QUANTILES[1])


@dataclass(frozen=True)
class Plan:
    """One (N, k): the program's resolved run and the options of its tasks."""

    run: glue.RRGRun
    options: tuple[str, ...]  # program arguments, without --tile-id


@dataclass(frozen=True)
class Array:
    """One sbatch array: array index ``i`` runs tile ``offset + i``."""

    jobname: str
    spec: str  # the --array value
    offset: int
    n_tasks: int
    script: str


# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
def check_options(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    """The serializer's own options: ranges, the values bound for
    ``#SBATCH`` lines, and options the chosen action would ignore."""
    _check_ranges(args)
    if args.action in SCRIPT_ACTIONS:
        check_sbatch_value("--array-time", args.array_time, required=False)
        check_sbatch_value(
            "--array-partition", args.array_partition, required=False
        )
        check_sbatch_value("the library path", str(LIB_ROOT))
    for dest, actions in OPTION_ACTIONS.items():
        given = getattr(args, dest) != parser.get_default(dest)
        if given and args.action not in actions:
            flag = f"--{dest.replace('_', '-')}"
            raise ValueError(f"{flag} does not apply to --{args.action}")
    if args.verify and not (args.action == "status" or args.only_missing):
        raise ValueError("--verify applies to --status and --only-missing")


def _check_ranges(args: argparse.Namespace) -> None:
    if min(args.N_list) < 1 or min(args.k_list) < 1:
        raise ValueError("--N-list and --k-list take positive integers")
    if args.array_max < 1 or args.array_concurrent < 1:
        raise ValueError("--array-max and --array-concurrent must be >= 1")
    if args.mem_mb < 0 or args.boot < 0:
        raise ValueError("--mem-mb and --boot must be >= 0")
    if not 0 <= args.boot_seed < TC_BOOT_SEED_LIMIT:
        raise ValueError(f"--boot-seed must lie in [0, {TC_BOOT_SEED_LIMIT})")
    window = args.peak_window
    if window and not window[0] < window[1]:  # also refuses a NaN bound
        raise ValueError("--peak-window needs LO < HI")


def check_forwarded(extras: Sequence[str]) -> None:
    """Refuse a forwarded option that neither program knows, naming the
    closest known one (a mistyped serializer option lands here). The
    program's flags come from its parser, so ``--no-verbose`` and the other
    generated forms count."""
    program = set(program_parser()._option_string_actions)
    ours = _flags(RRG_TransCluster_srun_args, RRG_TransCluster_srun_actions)
    for word in extras:
        flag = word.split("=", 1)[0]
        if not flag.startswith("-") or flag in program or _is_number(flag):
            continue
        close = difflib.get_close_matches(flag, sorted(program | ours), n=1)
        hint = f" (did you mean {close[0]}?)" if close else ""
        raise ValueError(f"unrecognized option {flag}{hint}")


def _flags(*specs: Mapping[tuple[str, ...], dict]) -> set[str]:
    return {flag for spec in specs for flags in spec for flag in flags}


def _is_number(word: str) -> bool:
    try:
        float(word)
    except ValueError:
        return False
    return True


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------
def program_parser() -> argparse.ArgumentParser:
    """The program's parser, named after it in its error messages."""
    parser = build_program_parser()
    parser.prog = PROGRAM.name
    return parser


def task_options(
    n_nodes: int, k: int, extras: Sequence[str]
) -> tuple[list[str], argparse.Namespace]:
    """The program arguments of (N, k) as typed, and the typed values of the
    ``PINNED`` options (``None`` if not typed); ``--tile-id`` is ours."""
    options = [str(n_nodes), str(k), *extras]
    parser = program_parser()
    parser.set_defaults(**dict.fromkeys((*PINNED, "tile_id")))
    typed = parser.parse_args(options)
    if typed.tile_id is not None:
        raise ValueError(f"{TILE_FLAG} is set per task by the serializer")
    return options, typed


def plan_run(n_nodes: int, k: int, extras: Sequence[str]) -> Plan:
    """Resolve (N, k) as the program does, then pin what was not typed."""
    options, typed = task_options(n_nodes, k, extras)
    args = program_parser().parse_args(options)
    if typed.ctrl_chunk is None and n_nodes >= TC_SPLIT_CTRL_MIN_N:
        args.ctrl_chunk = TC_SPLIT_CTRL_CHUNK
    run = glue.resolve_run(args, [str(PROGRAM), *options])
    resolved = {**vars(args), "sampler": run.method}
    for dest in PINNED:
        if resolved[dest] != getattr(typed, dest):
            options += [PROGRAM_FLAGS[dest], str(resolved[dest])]
    run = dataclasses.replace(run, argv=(str(PROGRAM), *options))
    return Plan(run, tuple(options))


def plan_runs(
    n_list: Sequence[int], k_list: Sequence[int], extras: Sequence[str]
) -> list[Plan]:
    """Resolve every (N, k) with the program's parser and ``resolve_run``."""
    check_forwarded(extras)
    plans = []
    for n_nodes in dict.fromkeys(n_list):
        for k in dict.fromkeys(k_list):
            try:
                plans.append(plan_run(n_nodes, k, extras))
            except ValueError as err:
                raise ValueError(f"N={n_nodes} k={k}: {err}") from err
    return plans


def tile_states(run: glue.RRGRun, verify: bool = False) -> dict[str, list[int]]:
    """Tile ids of ``run`` by state, judged by file name (``TILE_STATES``).
    ``verify`` reads every complete file as a task would, and also reports
    the ``unreadable`` ones (recomputed) and the ``foreign`` ones (the task
    exits with ``TC_EXIT_CFG_MISMATCH``)."""
    tiles_dir = run.rundir / TC_TILES_DIR
    names = set(os.listdir(tiles_dir)) if tiles_dir.is_dir() else set()
    header, ghash = glue.tile_header(run), grid_hash(run.grid)
    kinds = VERIFIED_STATES if verify else TILE_STATES
    states: dict[str, list[int]] = {name: [] for name in kinds}
    for tile_id in run.tile_ids:
        tile = decode_tile(run.layout, tile_id)
        paths = tile_paths(run.rundir, ghash, tile)
        if paths.final.name not in names:
            state = "partial" if paths.part.name in names else "missing"
        elif verify:
            state = _final_state(paths.final, header, tile)
        else:
            state = "complete"
        states[state].append(tile_id)
    return states


def _final_state(final: Path, header: TileHeader, tile: Tile) -> str:
    """What a task finds in a complete tile file."""
    try:
        check_tile(final, read_tile(final), header, tile, complete=True)
    except FileNotFoundError:  # set aside meanwhile
        return "missing"
    except TileUnreadable:
        return "unreadable"
    except TileMismatch:
        return "foreign"
    return "complete"


def pending_tiles(plan: Plan, args: argparse.Namespace) -> list[int]:
    """Every tile, or with ``--only-missing`` the tiles a task would compute
    (``RECOMPUTED``). Tiles of another run are reported, not resubmitted."""
    if not args.only_missing:
        return list(plan.run.tile_ids)
    states = tile_states(plan.run, args.verify)
    pending = sorted(i for state in RECOMPUTED for i in states.get(state, ()))
    if states.get("foreign"):
        print(
            f"{_label(plan)}: {len(states['foreign'])} tiles of another run "
            "are not resubmitted (see --status --verify)",
            file=sys.stderr,
        )
    elif not pending:
        print(f"{_label(plan)}: all tiles complete", file=sys.stderr)
    return pending


def task_mem_mb(n_nodes: int, engine: str, override: int) -> int:
    """``--mem`` of one task: the override, else the estimate in the defs."""
    if override:
        return override
    per_node = TC_MEM_MB_PER_NODE + TC_MEM_ENGINE_MB_PER_NODE[engine]
    return math.ceil(TC_MEM_SAFETY * (TC_MEM_BASE_MB + n_nodes * per_node))


def task_env() -> dict[str, str]:
    """Variables every task exports: thread caps, and the data root the run
    directories were resolved under."""
    threads = dict.fromkeys(SRUN_THREAD_ENV_VARS, str(DEFAULT_SRUN_THREADS))
    data = Path(glue.PATHDATA).absolute()  # the tasks run in LIB_ROOT
    return {**threads, **SBATCH_ENV, "LRGSG_DATA": str(data)}


def task_arrays(
    plan: Plan, args: argparse.Namespace, tile_ids: Sequence[int]
) -> list[Array]:
    """The sbatch arrays that run ``tile_ids`` of ``plan``."""
    run, arrays = plan.run, []
    command = [sys.executable, str(PROGRAM.relative_to(LIB_ROOT))]
    for offset, indices in split_arrays(tile_ids, args.array_max):
        jobname = build_jobname(
            program_short=RRG_TransCluster_progNameShrt,
            tokens=[
                f"N{run.n_nodes}",
                f"k{run.k}",
                f"a{offset // args.array_max}",
            ],
        )
        spec = format_array_spec(indices, args.array_concurrent)
        script = render_sbatch_tiles_script(
            jobname=jobname,
            chdir=str(LIB_ROOT),
            array=spec,
            mem_mb=task_mem_mb(run.n_nodes, run.engine, args.mem_mb),
            cpus=DEFAULT_SRUN_THREADS,
            time=args.array_time,
            partition=args.array_partition,
            signal=STOP_SIGNAL,
            env=task_env(),
            command=[*command, *plan.options],
            index_flag=TILE_FLAG,
            tile_offset=offset,
            log_dir=SBATCH_LOG_DIR,
        )
        arrays.append(Array(jobname, spec, offset, len(indices), script))
    return arrays


def _label(plan: Plan) -> str:
    return f"N={plan.run.n_nodes} k={plan.run.k}"


# ---------------------------------------------------------------------------
# Actions: (plans, args, command) -> exit code
# ---------------------------------------------------------------------------
def print_scripts(
    plans: Sequence[Plan], args: argparse.Namespace, command: Sequence[str]
) -> int:
    """Print the sbatch scripts. Nothing is written (``--only-missing``
    reads the tiles directory)."""
    for plan in plans:
        for array in task_arrays(plan, args, pending_tiles(plan, args)):
            print(array.script)
    return TC_EXIT_OK


def submit_arrays(
    plans: Sequence[Plan], args: argparse.Namespace, command: Sequence[str]
) -> int:
    """Publish and check every cfg.json and submissions/ directory, then
    submit the arrays and log each one. Nothing is submitted if a run
    directory holds another run or cannot take a log; after a failed sbatch
    or log, the arrays left out are listed."""
    if shutil.which("sbatch") is None:
        print("error: sbatch not found", file=sys.stderr)
        return TC_EXIT_FAILED
    for plan in plans:
        try:
            publish_cfg_sidecar(plan.run.rundir, **plan.run.identity)
        except FileExistsError as err:
            print(f"error: {err}", file=sys.stderr)
            return TC_EXIT_CFG_MISMATCH
    _warn_sbatch_env()
    (LIB_ROOT / SBATCH_LOG_DIR).mkdir(exist_ok=True)
    queue = [
        (plan, array)
        for plan in plans
        for array in task_arrays(plan, args, pending_tiles(plan, args))
    ]
    try:
        _check_log_dirs(queue)
    except OSError as err:
        print(f"error: {err}", file=sys.stderr)
        return TC_EXIT_FAILED
    return _submit_queue(queue, command)


def _check_log_dirs(queue: Sequence[tuple[Plan, Array]]) -> None:
    """Create ``submissions/`` in every run directory of ``queue`` and write
    a probe file there: a log that cannot be written raises ``OSError``
    before the first sbatch, not after it."""
    for rundir in dict.fromkeys(plan.run.rundir for plan, _ in queue):
        logs = rundir / TC_SUBMISSIONS_DIR
        logs.mkdir(exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=logs):
            pass


def _submit_queue(
    queue: Sequence[tuple[Plan, Array]], command: Sequence[str]
) -> int:
    """sbatch and log each array; its job line is printed as soon as sbatch
    returns. The first failure stops the queue and lists the arrays left
    out."""
    for done, (plan, array) in enumerate(queue):
        try:
            job_id = sbatch(array.script)
        except subprocess.CalledProcessError as err:
            print(f"error: sbatch: {err.stderr.strip()}", file=sys.stderr)
            _report_unsubmitted(queue[done:])
            return TC_EXIT_FAILED
        print(
            f"{array.jobname}: job {job_id}, {array.n_tasks} tasks", flush=True
        )
        try:
            log_submission(plan, array, job_id, command)
        except OSError as err:  # the next logs would fail alike: stop
            print(
                f"error: job {job_id} is queued, its log was not written: "
                f"{err}",
                file=sys.stderr,
            )
            if queue[done + 1 :]:
                _report_unsubmitted(queue[done + 1 :])
            return TC_EXIT_FAILED
    if queue:
        print(f"task logs: {LIB_ROOT / SBATCH_LOG_DIR}", file=sys.stderr)
    return TC_EXIT_OK


def _warn_sbatch_env() -> None:
    """sbatch lets its input variables override the #SBATCH lines."""
    names = sorted(n for n in os.environ if n.startswith(SBATCH_ENV_PREFIX))
    if names:
        print(
            f"warning: {', '.join(names)} set; sbatch lets them override "
            "the #SBATCH lines of the scripts",
            file=sys.stderr,
        )


def _report_unsubmitted(queue: Sequence[tuple[Plan, Array]]) -> None:
    print(
        "not submitted (--only-missing resubmits them once the queued "
        "arrays have ended):",
        file=sys.stderr,
    )
    for plan, array in queue:
        tiles = f"tiles {array.offset} + [{array.spec}]"
        print(f"  {_label(plan)} {array.jobname}: {tiles}", file=sys.stderr)


def sbatch(script: str) -> str:
    """Submit ``script``; return the job id."""
    result = subprocess.run(
        ["sbatch", "--parsable"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout.strip().split(";")[0]  # "jobid[;cluster]"


def log_submission(
    plan: Plan, array: Array, job_id: str, command: Sequence[str]
) -> Path:
    """Write ``submissions/<job id>.json`` in the run directory: the array,
    the layout its tile ids refer to, and the code it was submitted from."""
    run = plan.run
    entry = {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "job_id": job_id,
        "jobname": array.jobname,
        "array": array.spec,
        "offset": array.offset,
        "n_tasks": array.n_tasks,
        "grid": grid_hash(run.grid),
        "layout": dataclasses.asdict(run.layout),
        "rundir": str(run.rundir),
        "git": git_state(),
        "options": list(plan.options),
        "command": list(command),
    }
    data = json.dumps(entry, indent=1, sort_keys=True).encode()
    path = run.rundir / TC_SUBMISSIONS_DIR / f"{job_id}.json"
    return atomic_write_bytes(path, lambda fh: fh.write(data))


def report_status(
    plans: Sequence[Plan], args: argparse.Namespace, command: Sequence[str]
) -> int:
    """Tile counts per (N, k) of the grid and layout given."""
    for plan in plans:
        states = tile_states(plan.run, args.verify)
        counts = ", ".join(f"{len(ids)} {name}" for name, ids in states.items())
        print(f"{_label(plan)}: {counts} ({plan.run.rundir})")
    return TC_EXIT_OK


def collect_runs(
    plans: Sequence[Plan], args: argparse.Namespace, command: Sequence[str]
) -> int:
    """Aggregate each run directory; print its cost table and chi peak. A
    run that cannot be collected is reported and the next one tried."""
    code = TC_EXIT_OK
    for plan in plans:
        try:
            agg = collect(
                plan.run.rundir,
                plan.run.n_nodes,
                window=tuple(args.peak_window) if args.peak_window else None,
                n_boot=args.boot,
                seed=args.boot_seed,
                argv=command,
                legacy=args.legacy_txt,
            )
        except (TileError, ValueError, OSError) as err:
            print(f"{_label(plan)}: {err}", file=sys.stderr)
            code = TC_EXIT_FAILED
            continue
        print(f"{_label(plan)}: {plan.run.rundir / TC_AGG_NAME}")
        print(cost_table(agg["summary"]))
        meta = agg["meta"]
        _print_peak(meta)
        if meta["peak_at_edge"] == -1:  # agg.npz is still valid: exit 0
            reason = _no_peak_reason(agg["summary"], meta)
            print(f"{_label(plan)}: no chi peak: {reason}", file=sys.stderr)
        elif math.isnan(meta["peak_lo"]) and meta["n_boot"] > 0:
            reason = _no_interval_reason(meta)
            print(
                f"{_label(plan)}: no peak interval: {reason}", file=sys.stderr
            )
    return code


def _print_peak(meta) -> None:
    """The chi peak and its bootstrap interval, whose level the
    ``TC_BOOT_QUANTILES`` set. A maximum at a window end only bounds the
    location; ``EDGE_FRAC_WARN`` of the usable resamples peaking there may
    clip the interval."""
    level = TC_BOOT_QUANTILES[1] - TC_BOOT_QUANTILES[0]
    print(
        f"chi peak {meta['peak']:.4g} "
        f"[{meta['peak_lo']:.4g}, {meta['peak_hi']:.4g}] ({level:.0%}) "
        f"from {meta['peak_reps']} paired reps"
    )
    if meta["peak_at_edge"] == 1:
        warning = (
            "the chi maximum is at a window end; the location is only a bound"
        )
    elif meta["peak_edge_frac"] >= EDGE_FRAC_WARN:  # False if NaN: no resample
        warning = (
            f"{meta['peak_edge_frac']:.0%} of the usable resamples peak at a "
            "window end; the interval may be clipped"
        )
    else:
        return
    print(
        f"warning: {warning} (widen --peak-window or the grid)", file=sys.stderr
    )


def _no_peak_reason(summary, meta) -> str:
    """Why ``bootstrap_peak`` found no chi peak (flag -1): too few control
    values with converged reps in the window, or a chi that is flat or,
    below two paired reps, undefined."""
    ctrl, lo, hi = summary["ctrl"], meta["window_lo"], meta["window_hi"]
    n_ctrl = int(((ctrl >= lo) & (ctrl <= hi) & (summary["n_rep"] > 0)).sum())
    if n_ctrl < TC_PEAK_MIN_CTRLS:
        return (
            f"{n_ctrl} control values with converged reps in the window "
            f"(need {TC_PEAK_MIN_CTRLS})"
        )
    return (
        "chi is flat or undefined over the window "
        f"({meta['peak_reps']} paired reps)"
    )


def _no_interval_reason(meta) -> str:
    """Why a chi peak has no bootstrap interval although ``--boot`` > 0."""
    if meta["peak_reps"] < TC_BOOT_MIN_REPS:
        return f"{meta['peak_reps']} paired reps (need {TC_BOOT_MIN_REPS})"
    return "chi is flat in every resample"


ACTIONS = {
    "print": print_scripts,
    "exec": submit_arrays,
    "status": report_status,
    "collect": collect_runs,
}
