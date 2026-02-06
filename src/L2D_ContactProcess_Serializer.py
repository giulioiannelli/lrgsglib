from __future__ import annotations

import subprocess
from pathlib import Path

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    resolve_slanzarv_mode,
    _collect_values,
    _collect_values_typed,
    _determine_precision,
    _format_value_consistently,
    build_slanzarv_command,
    format_slanzarv_command,
)
from parsers.L2D_ContactProcess_Serializer import (
    L2D_CPROC_progname,
    L2D_CPROC_progname_shrt,
    parser,
)
from lrgsglib.config.progargs import (
    DEFAULT_L2D_CONTACTPROCESS_SRUN_GAMMA_LIST,
    DEFAULT_L2D_CONTACTPROCESS_SRUN_L_LIST,
    DEFAULT_L2D_CONTACTPROCESS_SRUN_MU_LIST,
    DEFAULT_L2D_CONTACTPROCESS_SRUN_P_LIST,
)


def main() -> None:
    args, unknown = parser.parse_known_args()
    exec_bool, print_bool = args.exec, args.print
    if not (exec_bool or print_bool):
        return

    use_slanzarv, dynamics_mode = resolve_slanzarv_mode(
        getattr(args, "mode", None),
        default_use_slanzarv=True,
    )

    dynamics = (dynamics_mode or "EI").upper()
    runlang = getattr(args, "runlang", "").upper() or "C1C"

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

    p_precision = _determine_precision(p_values)
    gamma_precision = _determine_precision(gamma_values)
    mu_precision = _determine_precision(mu_values)

    memoryfunc = build_memory_function(args.slanzarv_minMB, args.slanzarv_maxMB, L_values)
    script_path = Path("src") / f"{L2D_CPROC_progname}.py"

    total_printed = total_executed = 0

    def dispatch(L: int, pflip: float, gamma: float, mu: float) -> None:
        nonlocal total_executed, total_printed

        prog_args = [
            str(L),
            _format_value_consistently(pflip, p_precision),
            "--dynamics",
            dynamics,
            "--runlang",
            runlang,
        ]

        if dynamics == "EI":
            prog_args.extend([
                "--gamma",
                _format_value_consistently(gamma, gamma_precision),
            ])
        elif dynamics == "SIR":
            prog_args.extend([
                "--mu",
                _format_value_consistently(mu, mu_precision),
            ])

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
                f"p{_format_value_consistently(pflip, p_precision)}",
                f"dyn{dynamics}",
            ]
            if dynamics == "EI":
                jobname_tokens.append(f"ga{_format_value_consistently(gamma, gamma_precision)}")
            elif dynamics == "SIR":
                jobname_tokens.append(f"mu{_format_value_consistently(mu, mu_precision)}")

            jobname = build_jobname(
                program_short=L2D_CPROC_progname_shrt,
                tokens=jobname_tokens,
                job_id=args.slanzarv_id or None,
                prefix_job_id=True,
            )
            slanz_opts.extend(["--jobname", jobname])

            final_cmd = build_slanzarv_command(slanz_opts, cmd)

        if print_bool:
            if use_slanzarv:
                print(format_slanzarv_command(slanz_opts, cmd))
            else:
                print(" ".join(final_cmd))
            total_printed += 1
        if exec_bool:
            subprocess.run(final_cmd, check=True)
            total_executed += 1

    if dynamics == "EI":
        mu_seed = mu_values[0] if mu_values else 0.0
        for L in L_values:
            for pflip in p_values:
                for gamma in gamma_values:
                    dispatch(L, pflip, gamma, mu_seed)
    elif dynamics == "SIR":
        gamma_seed = gamma_values[0] if gamma_values else 0.0
        for L in L_values:
            for pflip in p_values:
                for mu in mu_values:
                    dispatch(L, pflip, gamma_seed, mu)
    else:
        raise ValueError(f"Unsupported dynamics selection for serializer: {dynamics}")

    if exec_bool:
        print(f"Total number of jobs executed: {total_executed}")
    if print_bool:
        print(f"Total number of jobs described: {total_printed}")


if __name__ == "__main__":
    main()
