"""Batch job serializer for L2D_VoterModel.py.

Sweeps lattice side x pflip x the voter dynamics axes (rule, update schedule,
nonlinearity alpha, noise eps, q-voter q) and submits one job per parameter
point via slanzarv, or one SLURM array per point (``--dispatch array``) with one
quench realisation per array task.

The voter dynamics axes default to a single value (the scalar program flags) and
become sweep axes only when their ``--*-list`` / ``--*-linsp`` flags are given.
Combinations are de-duplicated by their EFFECTIVE parameters, so e.g. sweeping
``--alpha-list`` while ``--rule linear`` does not emit duplicate jobs (alpha is
ignored by the linear rule).
"""

from __future__ import annotations

import subprocess
import sys
from itertools import product
from pathlib import Path

from tqdm import tqdm as _tqdm

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    build_slanzarv_command,
    format_slanzarv_command,
    render_sbatch_array_script,
)
from lrgsglib.config.progargs import L2D_VOTER_progname, L2D_VOTER_progname_shrt
from lrgsglib.statsys.VoterModel.defaults import (
    DEFAULT_ABSORBING_EVERY,
    DEFAULT_SOUT_EVERY,
    GILLESPIE_RULES,
    LINK_ONLY_RULES,
)
from parsers.L2D_VoterModel_Serializer import parser


def _axis(args, list_attr, scalar_attr, linsp_attr=None):
    """Resolve one sweep axis: explicit list, else linspace, else the scalar."""
    vals = getattr(args, list_attr, None)
    if vals:
        return list(vals)
    if linsp_attr is not None:
        linsp = getattr(args, linsp_attr, None)
        if linsp is not None and len(linsp):
            return list(linsp)
    return [getattr(args, scalar_attr)]


def _effective_key(rule, upd, alpha, eps, q):
    """Identify parameter combos that produce identical runs.

    Only the rule-relevant tunables enter the key: ``alpha``/``eps`` for the
    nonlinear rule, ``q``/``eps`` for the q-voter; linear and majority take none.
    """
    if rule == "nonlinear":
        return (rule, upd, round(float(alpha), 12), round(float(eps), 12))
    if rule == "qvoter":
        return (rule, upd, int(q), round(float(eps), 12))
    return (rule, upd)


def _resolve_combos(args):
    """Build the de-duplicated list of (rule, upd, alpha, eps, q) combos.

    Skips (and reports) combinations that are ill-defined: ``link`` and
    ``gillespie`` schedules require ``rule='linear'``.
    """
    rule_list = _axis(args, "rule_list", "rule")
    upd_list = _axis(args, "upd_mode_list", "upd_mode")
    alpha_list = _axis(args, "alpha_list", "alpha", "alpha_linsp")
    eps_list = _axis(args, "eps_list", "eps", "eps_linsp")
    q_list = _axis(args, "q_list", "q")

    seen, combos, skipped = set(), [], 0
    for rule, upd, alpha, eps, q in product(
        rule_list, upd_list, alpha_list, eps_list, q_list
    ):
        if upd == "link" and rule not in LINK_ONLY_RULES:
            skipped += 1
            continue
        if upd == "gillespie" and rule not in GILLESPIE_RULES:
            skipped += 1
            continue
        key = _effective_key(rule, upd, alpha, eps, q)
        if key in seen:
            continue
        seen.add(key)
        combos.append((rule, upd, alpha, eps, q))
    if skipped:
        print(
            f"note: skipped {skipped} (rule, upd_mode) combo(s) where "
            f"link/gillespie were paired with a non-linear rule.",
            file=sys.stderr,
        )
    return combos


def _combo_token(combo) -> str:
    """Compact, fully distinguishing jobname token for a combo."""
    rule, upd, alpha, eps, q = combo
    parts = [rule[:3]]
    if upd != "asynchronous":
        parts.append(
            {"synchronous": "sync", "link": "link", "gillespie": "gill"}.get(upd, upd)
        )
    if rule == "nonlinear":
        parts.append(f"a{alpha:.3g}")
        if eps:
            parts.append(f"e{eps:.3g}")
    elif rule == "qvoter":
        parts.append(f"q{int(q)}")
        if eps:
            parts.append(f"e{eps:.3g}")
    return "_".join(parts)


def _build_passthrough_flags(args, combo) -> list[str]:
    """Explicit CLI flags forwarded to L2D_VoterModel.py (no L / -p / -qid)."""
    rule, upd, alpha, eps, q = combo
    flags: list[str] = []

    # Graph engine / geometry / disorder cell
    ge = getattr(args, "graph_engine", "nx")
    if ge != "nx":
        flags += ["-ge", ge]
    flags += ["-g", str(args.geometry)]
    cell = getattr(args, "cell_type", "rand")
    if cell != "rand":
        flags += ["-c", cell]

    # Dynamics combo (Axis A rule + its tunables, Axis B schedule)
    flags += ["--rule", rule, "--upd_mode", upd]
    if rule == "nonlinear":
        flags += ["--alpha", f"{alpha:.6g}", "--eps", f"{eps:.6g}"]
    elif rule == "qvoter":
        flags += ["-q", str(int(q)), "--eps", f"{eps:.6g}"]

    # Backend / init / time
    flags += ["-rl", str(args.runlang), "-ic", str(args.init_cond),
              "-st", str(args.state_type)]
    if args.steps is not None:
        flags += ["--steps", str(args.steps)]
    if args.simref is not None:
        flags += ["--simref", f"{args.simref:.6g}"]
    # Number of realisations (array size in array mode; in-process loop otherwise)
    flags += ["-na", str(args.number_of_averages)]

    # Absorbing-state detection
    if args.absorbing_check is not None:
        flags += ["--absorbing-check"] if args.absorbing_check else ["--no-absorbing-check"]
    if args.absorbing_every != DEFAULT_ABSORBING_EVERY:
        flags += ["--absorbing-every", str(args.absorbing_every)]

    # Observable save gates
    if not args.save_magn:
        flags += ["--no-save-magn"]
    if args.save_sout:
        flags += ["--save-sout"]
    if args.save_cldist:
        flags += ["--save-cldist"]
    if args.save_all:
        flags += ["--save-all"]
    if args.save_cldist or args.save_all:
        flags += ["-cm", str(args.cluster_mode)]
    if args.save_sout or args.save_all:
        if args.sout_every != DEFAULT_SOUT_EVERY:
            flags += ["--sout-every", str(args.sout_every)]
        if args.sout_nlog is not None:
            flags += ["--sout-nlog", str(args.sout_nlog)]
        if args.sout_force_full:
            flags += ["--sout-force-full"]

    # Output naming / RNG
    if getattr(args, "out_suffix", ""):
        flags += ["-os", str(args.out_suffix)]
    if not args.randstr:
        flags += ["--no-randstr"]
    if args.seed is not None:
        flags += ["--seed", str(args.seed)]
    return flags


def main():
    args = parser.parse_args()
    exec_bool, print_bool = args.exec, getattr(args, "print", False)
    if not (exec_bool or print_bool):
        return

    progn = L2D_VOTER_progname
    progn_shrt = L2D_VOTER_progname_shrt

    side_list = list(args.side1_list)
    pflp_list = list(args.pflip_linsp)
    combos = _resolve_combos(args)

    memoryfunc = build_memory_function(
        args.slanzarv_minMB, args.slanzarv_maxMB, side_list
    )
    src_dir = Path(__file__).resolve().parent
    exec_path = src_dir.relative_to(Path.cwd()) / f"{progn}.py"

    total_printed = total_executed = 0
    dispatch_mode = getattr(args, "dispatch", "slanzarv")
    total_dispatch = len(side_list) * len(combos) * len(pflp_list)

    show_bar = bool(exec_bool) and total_dispatch > 1
    bar = (
        _tqdm(total=total_dispatch, desc=f"{dispatch_mode} submit", unit="job",
              file=sys.stderr, leave=True)
        if show_bar
        else None
    )

    def _advance_bar(jobname: str) -> None:
        if bar is not None:
            bar.set_postfix_str(jobname, refresh=False)
            bar.update(1)

    def dispatch_slanzarv(L, p, combo) -> None:
        nonlocal total_printed, total_executed
        passthrough = _build_passthrough_flags(args, combo)
        prog_args = [str(L), "-p", f"{p:.3g}"]
        cmd = ["python", str(exec_path), *prog_args, *passthrough]

        slanz_opts = ["-m", str(memoryfunc(L))]
        if getattr(args, "nomail", False):
            slanz_opts.append("--nomail")
        if getattr(args, "short", False):
            slanz_opts.append("--short")
        if getattr(args, "moretime", 0):
            slanz_opts.extend(["--time", str(args.moretime)])

        jobname = build_jobname(
            program_short=progn_shrt,
            tokens=[str(L), f"{p:.3g}", _combo_token(combo)],
            job_id=getattr(args, "slanzarv_id", None) or None,
            prefix_job_id=False,
        )
        slanz_opts.extend(["--jobname", jobname])
        slanz_cmd = build_slanzarv_command(slanz_opts, cmd)

        if print_bool:
            print(format_slanzarv_command(slanz_opts, cmd))
            total_printed += 1
        if exec_bool:
            subprocess.run(slanz_cmd)
            total_executed += 1
            _advance_bar(jobname)

    def dispatch_array(L, p, combo) -> None:
        nonlocal total_printed, total_executed
        passthrough = _build_passthrough_flags(args, combo)
        n_quench = int(args.number_of_averages)
        concurrent = int(min(getattr(args, "array_concurrent", n_quench), n_quench))
        jobname = build_jobname(
            program_short=progn_shrt,
            tokens=[str(L), f"{p:.3g}", _combo_token(combo), "arr"],
            job_id=getattr(args, "slanzarv_id", None) or None,
            prefix_job_id=False,
        )
        script = render_sbatch_array_script(
            jobname=jobname,
            array_range=(1, n_quench),
            concurrent=concurrent,
            mem_mb=int(memoryfunc(L)),
            partition=getattr(args, "array_partition", ""),
            time=getattr(args, "array_time", ""),
            exec_path=str(exec_path),
            L=int(L),
            p=float(p),
            passthrough=passthrough,
        )

        if print_bool:
            print(script)
            print()
            total_printed += 1
        if exec_bool:
            subprocess.run(["sbatch"], input=script, text=True, check=True)
            total_executed += 1
            _advance_bar(jobname)

    dispatch = dispatch_array if dispatch_mode == "array" else dispatch_slanzarv
    try:
        for L in side_list:
            for combo in combos:
                for p in pflp_list:
                    dispatch(L, p, combo)
    finally:
        if bar is not None:
            bar.close()

    print(f"Total number of jobs executed: {total_executed}")
    print(f"Total number of jobs printed: {total_printed}")


if __name__ == "__main__":
    main()
