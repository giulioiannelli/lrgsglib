"""Batch job serializer for L3D_IsingDynamics.py.

Generates slanzarv commands sweeping over lattice side lengths,
frustration fractions, and temperatures. All parameters are explicit
(no ``unknown`` passthrough).
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
from parsers.L3D_IsingDynamics import L3D_ISDYN_progname, L3D_ISDYN_progname_shrt
from parsers.L3D_IsingDynamics_Serialiser import parser


def _build_passthrough_flags(args) -> list[str]:
    """Build explicit CLI flags to forward to L3D_IsingDynamics.py."""
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

    # Ising dynamics
    rl = getattr(args, "runlang", "C1")
    flags.extend(["-rl", rl])

    # Averaging
    na = getattr(args, "number_of_averages", 1)
    flags.extend(["-na", str(na)])
    nt = getattr(args, "n_thermal", 1)
    if nt > 1:
        flags.extend(["--n_thermal", str(nt)])

    # Simulated Annealing
    if getattr(args, "sa_mode", False):
        flags.append("--sa-mode")
        flags.extend(["--T-init", str(args.T_init)])
        flags.extend(["--T-final", str(args.T_final)])
        flags.extend(["--cooling-schedule", args.cooling_schedule])
        flags.extend(["--cooling-rate", str(args.cooling_rate)])
        flags.extend(["--steps-per-T", str(args.steps_per_T)])
        flags.extend(["--n-temperatures", str(args.n_temperatures)])

    # Parallel Tempering
    if getattr(args, "pt_mode", False):
        flags.append("--pt-mode")
        flags.extend(["--n-replicas", str(args.n_replicas)])
        flags.extend(["--T-min", str(args.T_min)])
        flags.extend(["--T-max", str(args.T_max)])
        flags.extend(["--T-ladder-type", args.T_ladder_type])
        flags.extend(["--steps-per-exchange", str(args.steps_per_exchange)])
        flags.extend(["--n-exchanges", str(args.n_exchanges)])

    # Other Ising options
    ic = getattr(args, "init_cond", "ground_state_0")
    if ic != "ground_state_0":
        flags.extend(["-ic", ic])
    freq = getattr(args, "freq", 2)
    if freq != 2:
        flags.extend(["-fq", str(freq)])
    ts = getattr(args, "thrmsteps", 20)
    if ts != 20:
        flags.extend(["-ts", str(ts)])

    return flags


def main():
    args = parser.parse_args()
    exec_bool, print_bool = args.exec, getattr(args, "print", False)
    if not (exec_bool or print_bool):
        return

    progn = L3D_ISDYN_progname
    progn_shrt = L3D_ISDYN_progname_shrt

    side_list = list(args.side1_list)
    pflp_list = list(args.pflip_linsp)
    temp_list = list(args.Temp_linsp)

    memoryfunc = build_memory_function(
        args.slanzarv_minMB, args.slanzarv_maxMB, side_list
    )
    exec_path = LRGSG_SRC.relative_to(Path.cwd()) / f"{progn}.py"
    passthrough = _build_passthrough_flags(args)

    total_printed = total_executed = 0

    def dispatch(L: int, p: float, T: float) -> None:
        nonlocal total_printed, total_executed

        prog_args = [str(L), f"{p:.3g}", f"{T:.3g}"]
        cmd = ["python", str(exec_path), *prog_args, *passthrough]

        slanz_opts: list[str] = ["-m", str(memoryfunc(L))]
        if getattr(args, "nomail", False):
            slanz_opts.append("--nomail")
        if getattr(args, "short", False):
            slanz_opts.append("--short")
        if getattr(args, "moretime", 0):
            slanz_opts.extend(["--time", str(args.moretime)])

        jobname = build_jobname(
            program_short=progn_shrt,
            tokens=prog_args,
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

    for L in side_list:
        for p in pflp_list:
            for T in temp_list:
                dispatch(L, p, T)

    print(f"Total number of jobs executed: {total_executed}")
    print(f"Total number of jobs printed: {total_printed}")


if __name__ == "__main__":
    main()
