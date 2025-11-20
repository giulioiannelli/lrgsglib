from __future__ import annotations

import subprocess
from pathlib import Path

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    _collect_values,
    _collect_values_typed,
)
from parsers.L3D_TransCluster_Serialiser import parser, L3D_TransCluster_progName, L3D_TransCluster_progNameShrt
from lrgsglib.config.progargs import (
    DEFAULT_L3D_TRANSCLUSTER_SRUN_L_LIST,
    DEFAULT_L3D_TRANSCLUSTER_SRUN_P_LIST,
    DEFAULT_L3D_TRANSCLUSTER_SRUN_MU_LIST,
    DEFAULT_L3D_TRANSCLUSTER_SRUN_SIGMA_LIST,
)


def main() -> None:
    args, unknown = parser.parse_known_args()
    exec_bool, print_bool = args.exec, args.print
    if not (exec_bool or print_bool):
        return

    # Determine if we should use slanzarv based on mode
    use_slanzarv = False
    actual_mode = None
    if args.mode:
        if args.mode.startswith("slanzarv_"):
            use_slanzarv = True
            actual_mode = args.mode.replace("slanzarv_", "")
        else:
            actual_mode = args.mode

    L_values = _collect_values_typed(
        args.L_list,
        args.L_linsp,
        DEFAULT_L3D_TRANSCLUSTER_SRUN_L_LIST,
        int,
    )

    p_values = _collect_values(
        args.p_list,
        args.p_linsp,
        DEFAULT_L3D_TRANSCLUSTER_SRUN_P_LIST,
    )
    
    # Only collect mu/sigma values if explicitly provided (defaults are None)
    mu_values = None
    sigma_values = None
    if args.mu_list is not None or args.mu_linsp is not None:
        mu_values = _collect_values(
            args.mu_list,
            args.mu_linsp,
            DEFAULT_L3D_TRANSCLUSTER_SRUN_MU_LIST or [],
        )
    if args.sigma_list is not None or args.sigma_linsp is not None:
        sigma_values = _collect_values(
            args.sigma_list,
            args.sigma_linsp,
            DEFAULT_L3D_TRANSCLUSTER_SRUN_SIGMA_LIST or [],
        )

    memoryfunc = build_memory_function(args.slanzarv_minMB, args.slanzarv_maxMB, L_values)
    script_path = Path("src") / f"{L3D_TransCluster_progName}.py"

    # Determine dispatch mode:
    # - If p_list or p_linsp is explicitly provided, use "flip" mode (sweep over p)
    # - Otherwise, if mu or sigma lists are provided, use "normal" mode (sweep over mu/sigma)
    # - Default to flip mode if nothing is specified
    p_explicitly_provided = args.p_list is not None or args.p_linsp is not None
    use_normal_weights = (not p_explicitly_provided and 
                          (mu_values is not None or sigma_values is not None))

    total_printed = total_executed = 0

    def dispatch(L: int, probability: float, mu: float, sigma: float) -> None:
        nonlocal total_printed, total_executed

        # Build command: pass L and p as positional args,
        # and all other arguments via *unknown (passed through from command line)
        prog_args = [
            str(L),
            f"{probability:.12g}",
        ]

        # Add mode-specific flags
        if use_normal_weights:
            prog_args.extend([
                "--mu", f"{mu:.12g}",
                "--sigma", f"{sigma:.12g}",
                "--edge_weight", "normal"
            ])
        else:
            prog_args.extend(["--edge_weight", "flip"])

        # Add the actual mode if specified (without slanzarv_ prefix)
        if actual_mode:
            prog_args.extend(["--mode", actual_mode])

        cmd = ["python", str(script_path), *prog_args, *unknown]

        # Only build slanzarv command if use_slanzarv is True
        if use_slanzarv:
            slanz_opts = ["-m", str(memoryfunc(L))]
            if args.nomail:
                slanz_opts.append("--nomail")
            if args.short:
                slanz_opts.append("--short")
            if args.moretime:
                slanz_opts.extend(["--time", str(args.moretime)])

            # Build jobname from L, p, and optionally mu, sigma
            jobname_tokens = [f"L{L}", f"p{probability:.3g}"]
            if use_normal_weights:
                jobname_tokens.extend([f"m{mu:.3g}", f"s{sigma:.3g}"])
            
            jobname = build_jobname(
                program_short=L3D_TransCluster_progNameShrt,
                tokens=jobname_tokens,
                job_id=args.slanzarv_id or None,
                prefix_job_id=True,
            )
            slanz_opts.extend(["--jobname", jobname])

            final_cmd = ["slanzarv", *slanz_opts, *cmd]
        else:
            final_cmd = cmd

        if print_bool:
            print(" ".join(final_cmd))
            total_printed += 1
        if exec_bool:
            subprocess.run(final_cmd, check=True)
            total_executed += 1

    # Dispatch logic: if using normal weights, sweep over mu and sigma with fixed p;
    # otherwise sweep over p with fixed mu and sigma
    if use_normal_weights:
        p_seed = p_values[0] if p_values else 0.0
        # Use collected values or default to [0.0] if None
        mu_list = mu_values if mu_values is not None else [0.0]
        sigma_list = sigma_values if sigma_values is not None else [0.0]
        for L in L_values:
            for mu in mu_list:
                for sigma in sigma_list:
                    dispatch(L, p_seed, mu, sigma)
    else:
        mu_seed = 0.0
        sigma_seed = 0.0
        for L in L_values:
            for probability in p_values:
                dispatch(L, probability, mu_seed, sigma_seed)

    if exec_bool:
        print(f"Total number of jobs executed: {total_executed}")
    if print_bool:
        print(f"Total number of jobs described: {total_printed}")


if __name__ == "__main__":
    main()
