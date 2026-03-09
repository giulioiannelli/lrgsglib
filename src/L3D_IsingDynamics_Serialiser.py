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

    # Topological algorithm parameters
    rl_lower = rl.lower()
    if 'topo' in rl_lower:
        topo_n = getattr(args, 'topo_n_modes', 40)
        if topo_n != 40:
            flags.extend(['--topo-n-modes', str(topo_n)])
        sigma = getattr(args, 'topo_sigma_init', 0.15)
        if sigma != 0.15:
            flags.extend(['--topo-sigma-init', str(sigma)])
        chunk = getattr(args, 'topo_chunk_size', 5000)
        if chunk != 5000:
            flags.extend(['--topo-chunk-size', str(chunk)])
        if not getattr(args, 'topo_polish', True):
            flags.append('--no-topo-polish')
        psweeps = getattr(args, 'topo_polish_sweeps', 50)
        if psweeps != 50:
            flags.extend(['--topo-polish-sweeps', str(psweeps)])
        tau = getattr(args, 'topo_tau', 1.0)
        if tau != 1.0:
            flags.extend(['--topo-tau', str(tau)])
        fstr = getattr(args, 'topo_field_strength', 1.0)
        if fstr != 1.0:
            flags.extend(['--topo-field-strength', str(fstr)])
        # topo_fca needs SA params — ensure --sa-mode is set
        if 'topo_fca' in rl_lower and not getattr(args, 'sa_mode', False):
            flags.append('--sa-mode')
            flags.extend(['--T-init', str(getattr(args, 'T_init', 10.0))])
            flags.extend(['--T-final', str(getattr(args, 'T_final', 0.01))])
            flags.extend(['--cooling-schedule',
                          getattr(args, 'cooling_schedule', 'exponential')])
            flags.extend(['--cooling-rate',
                          str(getattr(args, 'cooling_rate', 0.95))])
            flags.extend(['--steps-per-T',
                          str(getattr(args, 'steps_per_T', 100))])
            flags.extend(['--n-temperatures',
                          str(getattr(args, 'n_temperatures', 100))])

    # Result persistence
    if getattr(args, 'save_results', False):
        flags.append('--save-results')
        sf = getattr(args, 'save_frequency', 0)
        if sf > 0:
            flags.extend(['-sf', str(sf)])
        # Observable selection (only forward non-defaults)
        if not getattr(args, 'save_ene', True):
            flags.append('--no-save-ene')
        if not getattr(args, 'save_magn', True):
            flags.append('--no-save-magn')
        if getattr(args, 'save_sout', False):
            flags.append('--save-sout')
        if getattr(args, 'save_all', False):
            flags.append('--save-all')

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
    es = getattr(args, "eqstep", 20)
    if es != 20:
        flags.extend(["-es", str(es)])

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

    # In SA/PT mode, temperature is controlled by the schedule — no T sweep
    use_sa = getattr(args, "sa_mode", False)
    use_pt = getattr(args, "pt_mode", False)
    if use_sa or use_pt:
        temp_list = [None]

    def dispatch(L: int, p: float, T: float | None) -> None:
        nonlocal total_printed, total_executed

        prog_args = [str(L), "-p", f"{p:.3g}"]
        jobname_tokens = [str(L), f"{p:.3g}"]
        if T is not None:
            prog_args.extend(["-T", f"{T:.3g}"])
            jobname_tokens.append(f"{T:.3g}")

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
            tokens=jobname_tokens,
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
