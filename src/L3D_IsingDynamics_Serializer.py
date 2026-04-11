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
from lrgsglib.config.progargs.defs.IsingDynamics import *
from parsers.L3D_IsingDynamics import L3D_ISDYN_progname, L3D_ISDYN_progname_shrt
from parsers.L3D_IsingDynamics_Serializer import parser


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
    cell = getattr(args, "cell_type", "rand")
    if cell != "rand":
        flags.extend(["-c", cell])
    pdil = getattr(args, "pdil", 0.0)
    if pdil:
        flags.extend(["--pdil", str(pdil)])
    ew = getattr(args, "edge_weight", "flip")
    if ew != "flip":
        flags.extend(["--edge_weight", ew])
        flags.extend(["--mu", str(getattr(args, "mu", 0.0))])
        flags.extend(["--sigma", str(getattr(args, "sigma", 0.0))])
    wd = getattr(args, "workdir", "")
    if wd:
        flags.extend(["-wd", wd])

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
        topo_n = getattr(args, 'topo_n_modes', DEFAULT_TOPO_N_MODES)
        if topo_n != DEFAULT_TOPO_N_MODES:
            flags.extend(['--topo-n-modes', str(topo_n)])
        sigma = getattr(args, 'topo_sigma_init', DEFAULT_TOPO_SIGMA_INIT)
        if sigma != DEFAULT_TOPO_SIGMA_INIT:
            flags.extend(['--topo-sigma-init', str(sigma)])
        chunk = getattr(args, 'topo_chunk_size', DEFAULT_TOPO_CHUNK_SIZE)
        if chunk != DEFAULT_TOPO_CHUNK_SIZE:
            flags.extend(['--topo-chunk-size', str(chunk)])
        if not getattr(args, 'topo_polish', DEFAULT_TOPO_POLISH):
            flags.append('--no-topo-polish')
        psweeps = getattr(args, 'topo_polish_sweeps', DEFAULT_TOPO_POLISH_SWEEPS)
        if psweeps != DEFAULT_TOPO_POLISH_SWEEPS:
            flags.extend(['--topo-polish-sweeps', str(psweeps)])
        tau = getattr(args, 'topo_tau', DEFAULT_TOPO_TAU)
        if tau != DEFAULT_TOPO_TAU:
            flags.extend(['--topo-tau', str(tau)])
        fstr = getattr(args, 'topo_field_strength', DEFAULT_TOPO_FIELD_STRENGTH)
        if fstr != DEFAULT_TOPO_FIELD_STRENGTH:
            flags.extend(['--topo-field-strength', str(fstr)])
        # topo_fca needs SA params — ensure --sa-mode is set
        if 'topo_fca' in rl_lower and not getattr(args, 'sa_mode', False):
            flags.append('--sa-mode')
            flags.extend(['--T-init', str(getattr(args, 'T_init', DEFAULT_SA_T_INIT))])
            flags.extend(['--T-final', str(getattr(args, 'T_final', DEFAULT_SA_T_FINAL))])
            flags.extend(['--cooling-schedule',
                          getattr(args, 'cooling_schedule', DEFAULT_SA_COOLING_SCHEDULE)])
            flags.extend(['--cooling-rate',
                          str(getattr(args, 'cooling_rate', DEFAULT_SA_COOLING_RATE))])
            flags.extend(['--steps-per-T',
                          str(getattr(args, 'steps_per_T', DEFAULT_SA_STEPS_PER_T))])
            flags.extend(['--n-temperatures',
                          str(getattr(args, 'n_temperatures', DEFAULT_SA_N_TEMPERATURES))])

    # CEM parameters
    if 'cem' in rl_lower:
        cem_iter = getattr(args, 'cem_iter', DEFAULT_CEM_ITER)
        if cem_iter != DEFAULT_CEM_ITER:
            flags.extend(['--cem-iter', str(cem_iter)])
        pop = getattr(args, 'cem_pop_size', DEFAULT_CEM_POP_SIZE)
        if pop != DEFAULT_CEM_POP_SIZE:
            flags.extend(['--cem-pop-size', str(pop)])
        elite = getattr(args, 'cem_elite_frac', DEFAULT_CEM_ELITE_FRAC)
        if elite != DEFAULT_CEM_ELITE_FRAC:
            flags.extend(['--cem-elite-frac', str(elite)])
        isig = getattr(args, 'cem_init_sigma', DEFAULT_CEM_INIT_SIGMA)
        if isig != DEFAULT_CEM_INIT_SIGMA:
            flags.extend(['--cem-init-sigma', str(isig)])
        sm = getattr(args, 'cem_smoothing', DEFAULT_CEM_SMOOTHING)
        if sm != DEFAULT_CEM_SMOOTHING:
            flags.extend(['--cem-smoothing', str(sm)])
        sf_ = getattr(args, 'cem_sigma_floor', DEFAULT_CEM_SIGMA_FLOOR)
        if sf_ != DEFAULT_CEM_SIGMA_FLOOR:
            flags.extend(['--cem-sigma-floor', str(sf_)])
        sc_ = getattr(args, 'cem_sigma_ceiling', DEFAULT_CEM_SIGMA_CEILING)
        if sc_ != DEFAULT_CEM_SIGMA_CEILING:
            flags.extend(['--cem-sigma-ceiling', str(sc_)])
        rst = getattr(args, 'cem_restarts', DEFAULT_CEM_RESTARTS)
        if rst != DEFAULT_CEM_RESTARTS:
            flags.extend(['--cem-restarts', str(rst)])
        if not getattr(args, 'cem_greedy', DEFAULT_CEM_GREEDY):
            flags.append('--no-cem-greedy')
        gsw = getattr(args, 'cem_greedy_sweeps', DEFAULT_CEM_GREEDY_SWEEPS)
        if gsw != DEFAULT_CEM_GREEDY_SWEEPS:
            flags.extend(['--cem-greedy-sweeps', str(gsw)])

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

    # Clustering
    nc = getattr(args, "NoClust", 1)
    if nc != 1:
        flags.extend(["-nc", str(nc)])
    vl = getattr(args, "val", "=1")
    if vl != "=1":
        flags.extend(["-vl", vl])

    # Output options
    os_ = getattr(args, "out_suffix", "")
    if os_:
        flags.extend(["-os", os_])
    if not getattr(args, "randstr", True):
        flags.append("--no-randstr")

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
    # Use __file__ to locate sibling program (robust across machines)
    src_dir = Path(__file__).resolve().parent
    exec_path = src_dir.relative_to(Path.cwd()) / f"{progn}.py"
    passthrough = _build_passthrough_flags(args)

    total_printed = total_executed = 0

    # In SA/PT/CEM mode, temperature is controlled internally — no T sweep
    use_sa = getattr(args, "sa_mode", False)
    use_pt = getattr(args, "pt_mode", False)
    use_cem = "cem" in getattr(args, "runlang", "").lower()
    if use_sa or use_pt or use_cem:
        temp_list = [None]

    def dispatch(L: int, p: float, T: float | None) -> None:
        nonlocal total_printed, total_executed

        prog_args = [str(L), "-p", f"{p:.3g}"]
        jobname_tokens = [str(L), f"{p:.3g}"]
        if T is not None:
            prog_args.extend(["-T", f"{T:.3g}"])
            jobname_tokens.append(f"{T:.3g}")
        if "topo" in getattr(args, "runlang", "").lower():
            jobname_tokens.append(f"M{getattr(args, 'topo_n_modes', DEFAULT_TOPO_N_MODES)}")

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
