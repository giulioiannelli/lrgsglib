"""Batch comparison: SA vs topo_met vs topo_fca on 3D lattices.

Generates 3x cluster jobs per (side, pflip) point:
  1. pb_sa          — standard simulated annealing (baseline)
  2. pb_topo_met    — topological Metropolis
  3. pb_topo_fca    — topological field cooling annealing

All save NPZ results via --save-results for post-analysis.

Usage
-----
Print commands without executing::

    python src/L3D_TopoComparison_Serializer.py -s1 8 10 -pFT 0.3 0.3 --print

Execute on cluster::

    python src/L3D_TopoComparison_Serializer.py -s1 8 10 12 \\
        -pFT 0.1 0.5 5 -na 100 --exec
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    build_slanzarv_command,
    format_slanzarv_command,
)
from lrgsglib import *
from lrgsglib.config.progargs.defs.IsingDynamics import DEFAULT_INIT_COND
from parsers.L3D_IsingDynamics import L3D_ISDYN_progname, L3D_ISDYN_progname_shrt
from parsers.L3D_TopoComparison_Serializer import parser

# Methods to compare: (runlang, needs_sa_mode, short_label)
_METHODS = [
    ("pb_sa",       True,  "SA"),
    ("pb_topo_met", False, "TM"),
    ("pb_topo_fca", False, "TF"),
]


def _build_method_flags(
    args, runlang: str, needs_sa_mode: bool
) -> list[str]:
    """Build CLI flags for a specific method."""
    flags: list[str] = []

    # Graph engine
    ge = getattr(args, "graph_engine", "nx")
    if ge != "nx":
        flags.extend(["-ge", ge])

    # Geometry
    geo = getattr(args, "geometry", "sc")
    if geo != "sc":
        flags.extend(["-g", geo])

    # Lattice parameters
    pdil = getattr(args, "pdil", 0.0)
    if pdil:
        flags.extend(["--pdil", str(pdil)])
    ew = getattr(args, "edge_weight", "flip")
    if ew != "flip":
        flags.extend(["--edge_weight", ew])
        flags.extend(["--mu", str(getattr(args, "mu", 0.0))])
        flags.extend(["--sigma", str(getattr(args, "sigma", 0.0))])

    # Runlang
    flags.extend(["-rl", runlang])

    # Averaging
    na = getattr(args, "number_of_averages", 100)
    flags.extend(["-na", str(na)])
    nt = getattr(args, "n_thermal", 1)
    if nt > 1:
        flags.extend(["--n_thermal", str(nt)])

    # SA parameters (needed for pb_sa and pb_topo_fca)
    if needs_sa_mode or "topo_fca" in runlang:
        flags.append("--sa-mode")
        flags.extend(["--T-init", str(getattr(args, "T_init", 10.0))])
        flags.extend(["--T-final", str(getattr(args, "T_final", 0.01))])
        flags.extend(["--cooling-schedule",
                       getattr(args, "cooling_schedule", "exponential")])
        flags.extend(["--cooling-rate",
                       str(getattr(args, "cooling_rate", 0.95))])
        flags.extend(["--steps-per-T",
                       str(getattr(args, "steps_per_T", 100))])
        flags.extend(["--n-temperatures",
                       str(getattr(args, "n_temperatures", 100))])

    # Topo parameters (for topo_met and topo_fca)
    if "topo" in runlang:
        topo_n = getattr(args, "topo_n_modes", 40)
        if topo_n != 40:
            flags.extend(["--topo-n-modes", str(topo_n)])
        sigma = getattr(args, "topo_sigma_init", 0.15)
        if sigma != 0.15:
            flags.extend(["--topo-sigma-init", str(sigma)])
        chunk = getattr(args, "topo_chunk_size", 5000)
        if chunk != 5000:
            flags.extend(["--topo-chunk-size", str(chunk)])
        if not getattr(args, "topo_polish", True):
            flags.append("--no-topo-polish")
        psweeps = getattr(args, "topo_polish_sweeps", 50)
        if psweeps != 50:
            flags.extend(["--topo-polish-sweeps", str(psweeps)])
        tau = getattr(args, "topo_tau", 1.0)
        if tau != 1.0:
            flags.extend(["--topo-tau", str(tau)])
        fstr = getattr(args, "topo_field_strength", 1.0)
        if fstr != 1.0:
            flags.extend(["--topo-field-strength", str(fstr)])

    # Other Ising options
    ic = getattr(args, "init_cond", DEFAULT_INIT_COND)
    if ic != DEFAULT_INIT_COND:
        flags.extend(["-ic", ic])
    freq = getattr(args, "freq", 2)
    if freq != 2:
        flags.extend(["-fq", str(freq)])
    ts = getattr(args, "thrmsteps", 20)
    if ts != 20:
        flags.extend(["-ts", str(ts)])

    # Always save results
    flags.append("--save-results")

    return flags


def main():
    args = parser.parse_args()
    exec_bool = args.exec
    print_bool = getattr(args, "print", False)
    if not (exec_bool or print_bool):
        return

    progn_shrt = L3D_ISDYN_progname_shrt

    side_list = list(args.side1_list)
    pflp_list = list(args.pflip_linsp)

    memoryfunc = build_memory_function(
        args.slanzarv_minMB, args.slanzarv_maxMB, side_list
    )
    src_dir = Path(__file__).resolve().parent
    exec_path = src_dir.relative_to(Path.cwd()) / f"{L3D_ISDYN_progname}.py"

    total_printed = total_executed = 0

    def dispatch(L: int, p: float, method_flags: list[str],
                 method_label: str) -> None:
        nonlocal total_printed, total_executed

        # For SA/TFCA the temperature is set by SA schedule, use T=0.5 as dummy
        prog_args = [str(L), "-p", f"{p:.3g}", "-T", "0.5"]
        cmd = ["python", str(exec_path), *prog_args, *method_flags]

        slanz_opts: list[str] = ["-m", str(memoryfunc(L))]
        if getattr(args, "nomail", False):
            slanz_opts.append("--nomail")
        if getattr(args, "short", False):
            slanz_opts.append("--short")
        if getattr(args, "moretime", 0):
            slanz_opts.extend(["--time", str(args.moretime)])

        jobname = build_jobname(
            program_short=progn_shrt,
            tokens=[str(L), f"{p:.3g}", method_label],
            job_id=getattr(args, "slanzarv_id", None) or None,
            prefix_job_id=False,
        )
        slanz_opts.extend(["--jobname", jobname])

        slanz_cmd = build_slanzarv_command(slanz_opts, cmd)

        if print_bool:
            print(format_slanzarv_command(slanz_opts, cmd))
            print()
            total_printed += 1
        if exec_bool:
            subprocess.run(slanz_cmd)
            total_executed += 1

    for L in side_list:
        for p in pflp_list:
            for runlang, needs_sa, label in _METHODS:
                flags = _build_method_flags(args, runlang, needs_sa)
                dispatch(L, p, flags, label)

    print(f"Total number of jobs executed: {total_executed}")
    print(f"Total number of jobs printed: {total_printed}")


if __name__ == "__main__":
    main()
