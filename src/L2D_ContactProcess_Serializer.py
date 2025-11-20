from __future__ import annotations

import subprocess
from pathlib import Path

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    resolve_slanzarv_mode,
    _determine_precision,
    _format_value_consistently,
    _collect_values,
    _collect_values_typed,
)
from parsers.L2D_ContactProcess_Serializer import parser, L2D_CPROC_progname, L2D_CPROC_progname_shrt
from lrgsglib.config.progargs import (
    DEFAULT_L2D_CONTACTPROCESS_SRUN_L_LIST,
    DEFAULT_L2D_CONTACTPROCESS_SRUN_P_LIST,
    DEFAULT_L2D_CONTACTPROCESS_SRUN_GAMMA_LIST,
    DEFAULT_L2D_CONTACTPROCESS_SRUN_MU_LIST,
    DEFAULT_DYNAMICS,
)


def main() -> None:
    args, unknown = parser.parse_known_args()
    exec_bool, print_bool = args.exec, args.print
    if not (exec_bool or print_bool):
        return

    use_slanzarv, _ = resolve_slanzarv_mode(getattr(args, "mode", None), default_use_slanzarv=True)

    L_values = _collect_values_typed(
        args.L_list,
        args.L_linsp,
        DEFAULT_L2D_CONTACTPROCESS_SRUN_L_LIST,
        int,
    )

    p_values = _collect_values(
        args.p_list,
        args.p_linsp,
        DEFAULT_L2D_CONTACTPROCESS_SRUN_P_LIST,
    )

    gamma_values = _collect_values(
        args.gamma_list,
        args.gamma_linsp,
        DEFAULT_L2D_CONTACTPROCESS_SRUN_GAMMA_LIST,
    )

    mu_values = _collect_values(
        args.mu_list,
        args.mu_linsp,
        DEFAULT_L2D_CONTACTPROCESS_SRUN_MU_LIST,
    )

    dynamics_values = [d.upper() for d in (args.dynamics_list or [DEFAULT_DYNAMICS])]

    p_precision = _determine_precision(p_values)
    gamma_precision = _determine_precision(gamma_values)
    mu_precision = _determine_precision(mu_values)

    memoryfunc = build_memory_function(args.slanzarv_minMB, args.slanzarv_maxMB, L_values)
    script_path = Path("src") / f"{L2D_CPROC_progname}.py"

    total_printed = total_executed = 0

    def dispatch(L: int, p_value: float, dynamics: str, rate_value: float) -> None:
        nonlocal total_printed, total_executed

        p_str = _format_value_consistently(p_value, p_precision)
        prog_args = [
            str(L),
            p_str,
            "--dynamics",
            dynamics,
        ]

        if dynamics == "EI":
            rate_str = _format_value_consistently(rate_value, gamma_precision)
            prog_args.extend(["--gamma", rate_str])
        elif dynamics == "SIR":
            rate_str = _format_value_consistently(rate_value, mu_precision)
            prog_args.extend(["--mu", rate_str])
        else:
            return

        cmd = ["python", str(script_path), *prog_args, *unknown]

        final_cmd = cmd
        if use_slanzarv:
            slanz_opts = ["-m", str(memoryfunc(L))]
            if args.nomail:
                slanz_opts.append("--nomail")
            if args.short:
                slanz_opts.append("--short")
            if args.moretime:
                slanz_opts.extend(["--time", str(args.moretime)])

            jobname_tokens = [
                f"L{L}",
                f"p{p_str}",
                dynamics,
                f"rate{rate_str}",
            ]
            jobname = build_jobname(
                program_short=L2D_CPROC_progname_shrt,
                tokens=jobname_tokens,
                job_id=args.slanzarv_id or None,
                prefix_job_id=True,
            )
            slanz_opts.extend(["--jobname", jobname])

            final_cmd = ["slanzarv", *slanz_opts, *cmd]

        if print_bool:
            print(" ".join(final_cmd))
            total_printed += 1
        if exec_bool:
            subprocess.run(final_cmd, check=True)
            total_executed += 1

    for L in L_values:
        for p_value in p_values:
            for dynamics in dynamics_values:
                if dynamics == "EI":
                    for gamma in gamma_values:
                        dispatch(L, p_value, dynamics, gamma)
                elif dynamics == "SIR":
                    for mu in mu_values:
                        dispatch(L, p_value, dynamics, mu)

    if exec_bool:
        print(f"Total number of jobs executed: {total_executed}")
    if print_bool:
        print(f"Total number of jobs described: {total_printed}")


if __name__ == "__main__":
    main()
