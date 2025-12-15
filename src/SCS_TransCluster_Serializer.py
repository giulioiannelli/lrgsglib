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
    build_slanzarv_command,
    format_slanzarv_command,
)
from parsers.SCS_TransCluster_Serializer import parser, SCS_TransCluster_progName, SCS_TransCluster_progNameShrt
from lrgsglib.config.progargs import (
    DEFAULT_SCS_TRANSCLUSTER_SRUN_N_LIST,
    DEFAULT_SCS_TRANSCLUSTER_SRUN_GAMMA_LIST,
    DEFAULT_SCS_TRANSCLUSTER_SRUN_J0_LIST,
    DEFAULT_SCS_TRANSCLUSTER_SRUN_J_LIST,
    DEFAULT_SCS_TRANSCLUSTER_SRUN_G_LIST,
)


def main() -> None:
    args, unknown = parser.parse_known_args()
    exec_bool, print_bool = args.exec, args.print
    if not (exec_bool or print_bool):
        return

    use_slanzarv, actual_mode = resolve_slanzarv_mode(
        getattr(args, "mode", None),
        default_use_slanzarv=True,
    )

    N_values = _collect_values_typed(
        args.N_list,
        None,  # No linspace for N
        DEFAULT_SCS_TRANSCLUSTER_SRUN_N_LIST,
        int,
    )

    gamma_values = _collect_values_typed(
        args.gamma_list,
        None,  # No linspace for gamma
        DEFAULT_SCS_TRANSCLUSTER_SRUN_GAMMA_LIST,
        float,
    )

    j0_values = _collect_values(
        args.J0_list,
        args.J0_linsp,
        DEFAULT_SCS_TRANSCLUSTER_SRUN_J0_LIST,
    )

    J_values = _collect_values(
        args.J_list,
        args.J_linsp,
        DEFAULT_SCS_TRANSCLUSTER_SRUN_J_LIST,
    )

    g_values = _collect_values(
        args.g_list,
        args.g_linsp,
        DEFAULT_SCS_TRANSCLUSTER_SRUN_G_LIST,
    )


    # Diagonal values are strings (use explicit values from argparse, which include defaults)
    diagonal_values: list[str] = [str(v) for v in args.diagonal_list]

    # Calculate precision needed for consistent formatting in jobnames
    gamma_precision = _determine_precision(gamma_values)
    j0_precision = _determine_precision(j0_values)
    J_precision = _determine_precision(J_values)
    g_precision = _determine_precision(g_values)

    memoryfunc = build_memory_function(args.slanzarv_minMB, args.slanzarv_maxMB, N_values)
    script_path = Path("src") / f"{SCS_TransCluster_progName}.py"

    total_printed = total_executed = 0

    def dispatch(N: int, gamma: float, j0: float, J_value: float, g_value: float, diagonal: str) -> None:
        nonlocal total_printed, total_executed

        # Build command: pass N and gamma as positional args,
        # and all other arguments via *unknown (passed through from command line)
        prog_args = [
            str(N),
            "--gamma",
            _format_value_consistently(gamma, gamma_precision),
            "--J0",
            _format_value_consistently(j0, j0_precision),
            "--J",
            _format_value_consistently(J_value, J_precision),
            "--g",
            _format_value_consistently(g_value, g_precision),
            "--diagonal",
            diagonal,
        ]

        if actual_mode:
            prog_args.extend(["--mode", actual_mode])

        cmd = ["python", str(script_path), *prog_args, *unknown]

        final_cmd = cmd
        if use_slanzarv:
            slanz_opts = ["-m", str(memoryfunc(N))]
            if args.nomail:
                slanz_opts.append("--nomail")
            if args.short:
                slanz_opts.append("--short")
            if args.moretime:
                slanz_opts.extend(["--time", str(args.moretime)])

            # Build jobname from N, gamma, J, g, j0, diagonal
            # Use consistent formatting based on calculated precision
            jobname_tokens = [
                f"N{N}",
                f"g{_format_value_consistently(gamma, gamma_precision)}",
                f"J{_format_value_consistently(J_value, J_precision)}",
                f"gcpl{_format_value_consistently(g_value, g_precision)}",
                f"J0{_format_value_consistently(j0, j0_precision)}",
                f"diag{diagonal}",
            ]
            jobname = build_jobname(
                program_short=SCS_TransCluster_progNameShrt,
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

    for N in N_values:
        for gamma in gamma_values:
            for j0 in j0_values:
                for J_value in J_values:
                    for g_value in g_values:
                        for diagonal in diagonal_values:
                            dispatch(N, gamma, j0, J_value, g_value, diagonal)

    if exec_bool:
        print(f"Total number of jobs executed: {total_executed}")
    if print_bool:
        print(f"Total number of jobs described: {total_printed}")


if __name__ == "__main__":
    main()
