"""Batch job serializer for L3D_IsingDynamics.py.

Generates slanzarv commands sweeping over lattice side lengths,
frustration fractions, and temperatures. All parameters are explicit
(no ``unknown`` passthrough).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from tqdm import tqdm as _tqdm

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    build_slanzarv_command,
    format_slanzarv_command,
    render_sbatch_array_script,
)
from lrgsglib import *
from lrgsglib.config.progargs.defs.IsingDynamics import *
from parsers.L3D_IsingDynamics import L3D_ISDYN_progname, L3D_ISDYN_progname_shrt
from parsers.L3D_IsingDynamics_Serializer import parser


def _apply_overrides(flags: list[str], overrides: dict[str, str] | None) -> list[str]:
    """Replace or inject ``--flag value`` pairs from an overrides mapping.

    Every key in ``overrides`` is a long flag string (e.g. ``"--topo-n-modes"``);
    if present in ``flags`` its adjacent value is replaced in place, otherwise
    the ``(flag, value)`` pair is appended.
    """
    if not overrides:
        return flags
    out = list(flags)
    for key, value in overrides.items():
        try:
            idx = out.index(key)
        except ValueError:
            out.extend([key, str(value)])
        else:
            if idx + 1 < len(out):
                out[idx + 1] = str(value)
            else:
                out.append(str(value))
    return out


def _build_passthrough_flags(
    args,
    *,
    overrides: dict[str, str] | None = None,
) -> list[str]:
    """Build explicit CLI flags to forward to L3D_IsingDynamics.py.

    ``overrides`` may supply replacements for per-sweep flags (e.g. one
    ``--topo-n-modes`` or ``--cem-iter`` value per array submission); it
    is applied after the args-derived flags are assembled.
    """
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
    ic = getattr(args, "init_cond", DEFAULT_INIT_COND)
    if ic != DEFAULT_INIT_COND:
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

    return _apply_overrides(flags, overrides)


def _resolve_sweep_lists(args) -> tuple[list[int | None], list[int | None]]:
    """Return ``(M_list, cem_iter_list)`` aligned pairwise.

    Falls back to the single ``--topo-n-modes`` / ``--cem-iter`` scalar
    when the corresponding list flag is absent. Entries are ``None`` when
    the underlying parameter should not override the scalar passthrough
    (i.e. single-value historical behaviour).
    """
    M_list_raw = getattr(args, "topo_n_modes_list", None)
    C_list_raw = getattr(args, "cem_iter_list", None)
    if M_list_raw is not None:
        M_list: list[int | None] = list(M_list_raw)
    else:
        M_list = [None]
    if C_list_raw is not None:
        if M_list_raw is None:
            raise SystemExit(
                "error: --cem-iter-list/-cL requires --topo-n-modes-list/-mL"
            )
        if len(C_list_raw) != len(M_list_raw):
            raise SystemExit(
                f"error: -cL length ({len(C_list_raw)}) must equal "
                f"-mL length ({len(M_list_raw)})"
            )
        iter_list: list[int | None] = list(C_list_raw)
    else:
        iter_list = [None] * len(M_list)
    return M_list, iter_list


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
    M_list, iter_list = _resolve_sweep_lists(args)

    memoryfunc = build_memory_function(
        args.slanzarv_minMB, args.slanzarv_maxMB, side_list
    )
    # Use __file__ to locate sibling program (robust across machines)
    src_dir = Path(__file__).resolve().parent
    exec_path = src_dir.relative_to(Path.cwd()) / f"{progn}.py"

    total_printed = total_executed = 0

    # In SA/PT/CEM mode, temperature is controlled internally — no T sweep
    use_sa = getattr(args, "sa_mode", False)
    use_pt = getattr(args, "pt_mode", False)
    use_cem = "cem" in getattr(args, "runlang", "").lower()
    if use_sa or use_pt or use_cem:
        temp_list = [None]

    dispatch_mode = getattr(args, "dispatch", DEFAULT_DISPATCH)

    # Pre-compute total dispatch count for the progress bar.
    if dispatch_mode == "array":
        total_dispatch = (
            len(side_list) * len(M_list) * len(pflp_list)
        )
    else:
        total_dispatch = (
            len(side_list) * len(M_list) * len(pflp_list) * len(temp_list)
        )

    # Progress bar: enabled whenever actually submitting (--exec). For
    # --print we stay silent so the heredoc stream is clean for piping.
    show_bar = bool(exec_bool) and total_dispatch > 1
    bar = (
        _tqdm(
            total=total_dispatch,
            desc=f"{dispatch_mode} submit",
            unit="job",
            file=sys.stderr,
            leave=True,
        )
        if show_bar
        else None
    )

    def _advance_bar(jobname: str) -> None:
        if bar is not None:
            bar.set_postfix_str(jobname, refresh=False)
            bar.update(1)

    def _overrides_for(M: int | None, n_iter: int | None) -> dict[str, str] | None:
        ov: dict[str, str] = {}
        if M is not None:
            ov["--topo-n-modes"] = str(M)
        if n_iter is not None:
            ov["--cem-iter"] = str(n_iter)
        return ov or None

    def _topo_mode_for_jobname(M: int | None) -> int:
        # Matches pre-existing behaviour: always display the active M in
        # the jobname for topo runlangs, defaulting to the scalar flag.
        if M is not None:
            return M
        return getattr(args, "topo_n_modes", DEFAULT_TOPO_N_MODES)

    def dispatch_slanzarv(
        L: int,
        p: float,
        T: float | None,
        M: int | None,
        n_iter: int | None,
    ) -> None:
        nonlocal total_printed, total_executed

        passthrough = _build_passthrough_flags(
            args, overrides=_overrides_for(M, n_iter),
        )
        prog_args = [str(L), "-p", f"{p:.3g}"]
        jobname_tokens = [str(L), f"{p:.3g}"]
        if T is not None:
            prog_args.extend(["-T", f"{T:.3g}"])
            jobname_tokens.append(f"{T:.3g}")
        if "topo" in getattr(args, "runlang", "").lower():
            jobname_tokens.append(f"M{_topo_mode_for_jobname(M)}")

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
            _advance_bar(jobname)

    def dispatch_array(
        L: int,
        p: float,
        M: int | None,
        n_iter: int | None,
    ) -> None:
        nonlocal total_printed, total_executed

        passthrough = _build_passthrough_flags(
            args, overrides=_overrides_for(M, n_iter),
        )
        n_quench = int(getattr(args, "number_of_averages", 1))
        concurrent = int(min(getattr(args, "array_concurrent", DEFAULT_ARRAY_CONCURRENT),
                             n_quench))
        jobname_tokens = [str(L), f"{p:.3g}"]
        if "topo" in getattr(args, "runlang", "").lower():
            jobname_tokens.append(f"M{_topo_mode_for_jobname(M)}")
        if n_iter is not None:
            jobname_tokens.append(f"I{n_iter}")
        jobname_tokens.append("arr")
        jobname = build_jobname(
            program_short=progn_shrt,
            tokens=jobname_tokens,
            job_id=getattr(args, "slanzarv_id", None) or None,
            prefix_job_id=False,
        )
        script = render_sbatch_array_script(
            jobname=jobname,
            array_range=(1, n_quench),
            concurrent=concurrent,
            mem_mb=int(memoryfunc(L)),
            partition=getattr(args, "array_partition", DEFAULT_ARRAY_PARTITION),
            time=getattr(args, "array_time", DEFAULT_ARRAY_TIME),
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

    try:
        if dispatch_mode == "array":
            for L in side_list:
                for M, n_iter in zip(M_list, iter_list):
                    for p in pflp_list:
                        dispatch_array(L, p, M, n_iter)
        else:
            for L in side_list:
                for M, n_iter in zip(M_list, iter_list):
                    for p in pflp_list:
                        for T in temp_list:
                            dispatch_slanzarv(L, p, T, M, n_iter)
    finally:
        if bar is not None:
            bar.close()

    print(f"Total number of jobs executed: {total_executed}")
    print(f"Total number of jobs printed: {total_printed}")


if __name__ == "__main__":
    main()
