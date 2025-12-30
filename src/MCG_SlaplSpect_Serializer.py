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
from parsers.MCG_SlaplSpect_Serializer import (
    MCG_MCGSS_progname,
    MCG_MCGSS_progname_shrt,
    parser,
)
from lrgsglib.config.progargs import (
    DEFAULT_MCG_SLAPLSPECT_SRUN_MATRIX_LIST,
    DEFAULT_MCG_SLAPLSPECT_SRUN_FRACTION_LIST,
    DEFAULT_MCG_SLAPLSPECT_SRUN_ITERATIONS_LIST,
    DEFAULT_MCG_SLAPLSPECT_SRUN_P_LIST,
)


def parse_matrix_list(matrix_str: str) -> list[tuple[float, float, float, float]]:
    """
    Parse matrix list from string format.

    Format: "p1,p2,p3,p4;p1,p2,p3,p4;..."
    Example: "0.8,0.6,0.6,0.8;0.9,0.5,0.5,0.9"

    Returns list of (p1, p2, p3, p4) tuples.
    """
    matrices = []
    for matrix_part in matrix_str.split(';'):
        values = [float(x.strip()) for x in matrix_part.split(',')]
        if len(values) != 4:
            raise ValueError(f"Each matrix must have exactly 4 values, got {len(values)}")
        matrices.append(tuple(values))
    return matrices


def parse_matrix_file(filepath: str) -> list[tuple[float, float, float, float]]:
    """
    Parse matrix list from file.

    File format (one matrix per line):
    p1 p2 p3 p4
    p1 p2 p3 p4
    ...

    Example file content:
    0.8 0.6 0.6 0.8
    0.9 0.5 0.5 0.9
    0.7 0.7 0.7 0.7

    Returns list of (p1, p2, p3, p4) tuples.
    """
    matrices = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):  # Skip empty lines and comments
                continue
            values = [float(x) for x in line.split()]
            if len(values) != 4:
                raise ValueError(f"Each line must have exactly 4 values, got {len(values)}")
            matrices.append(tuple(values))
    return matrices


def main() -> None:
    args, unknown = parser.parse_known_args()
    exec_bool, print_bool = args.exec, args.print

    if not (exec_bool or print_bool):
        return

    use_slanzarv, _ = resolve_slanzarv_mode(
        getattr(args, "mode", None),
        default_use_slanzarv=True,
    )

    # Parse cascade probability matrices
    if args.matrix_list is not None:
        matrices = parse_matrix_list(args.matrix_list)
    elif args.matrix_file is not None:
        matrices = parse_matrix_file(args.matrix_file)
    else:
        matrices = DEFAULT_MCG_SLAPLSPECT_SRUN_MATRIX_LIST

    # Collect other parameter sweep values
    fraction_values = _collect_values(
        args.fraction_list,
        args.fraction_linsp,
        DEFAULT_MCG_SLAPLSPECT_SRUN_FRACTION_LIST,
    )
    iterations_values = _collect_values_typed(
        args.iterations_list,
        args.iterations_linsp,
        DEFAULT_MCG_SLAPLSPECT_SRUN_ITERATIONS_LIST,
        int,
    )
    pflip_values = _collect_values(
        args.p_list,
        args.p_linsp,
        DEFAULT_MCG_SLAPLSPECT_SRUN_P_LIST,
    )

    # Determine precision for consistent formatting
    # Collect all p-values for precision determination
    all_p_values = [p for matrix in matrices for p in matrix]
    p_precision = _determine_precision(all_p_values)
    fraction_precision = _determine_precision(fraction_values)
    pflip_precision = _determine_precision(pflip_values)

    # Memory function (based on iterations, which controls graph size)
    memoryfunc = build_memory_function(
        args.slanzarv_minMB,
        args.slanzarv_maxMB,
        iterations_values
    )

    script_path = Path("src") / f"{MCG_MCGSS_progname}.py"

    total_printed = total_executed = 0

    def dispatch(p1: float, p2: float, p3: float, p4: float,
                 fraction: float, iterations: int, pflip: float) -> None:
        """Submit/print a single job with given parameters."""
        nonlocal total_executed, total_printed

        # Build program arguments
        prog_args = [
            _format_value_consistently(p1, p_precision),
            _format_value_consistently(p2, p_precision),
            _format_value_consistently(p3, p_precision),
            _format_value_consistently(p4, p_precision),
            _format_value_consistently(fraction, fraction_precision),
            str(iterations),
            "-p",
            _format_value_consistently(pflip, pflip_precision),
        ]

        cmd = ["python", str(script_path), *prog_args, *unknown]

        final_cmd = cmd
        if use_slanzarv:
            slanz_opts = ["-m", str(memoryfunc(iterations))]
            if args.nomail:
                slanz_opts.append("--nomail")
            if args.short:
                slanz_opts.append("--short")
            if args.moretime:
                slanz_opts.extend(["--time", str(args.moretime)])

            # Build job name with matrix hash for uniqueness
            matrix_hash = f"{_format_value_consistently(p1, p_precision)}" \
                         f"{_format_value_consistently(p2, p_precision)}" \
                         f"{_format_value_consistently(p3, p_precision)}" \
                         f"{_format_value_consistently(p4, p_precision)}"

            jobname_tokens = [
                f"M{matrix_hash}",
                f"fr{_format_value_consistently(fraction, fraction_precision)}",
                f"it{iterations}",
                f"pf{_format_value_consistently(pflip, pflip_precision)}",
            ]

            jobname = build_jobname(
                program_short=MCG_MCGSS_progname_shrt,
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

    # Parameter sweep: nested loops over all combinations
    for (p1, p2, p3, p4) in matrices:
        for fraction in fraction_values:
            for iterations in iterations_values:
                for pflip in pflip_values:
                    dispatch(p1, p2, p3, p4, fraction, iterations, pflip)

    if exec_bool:
        print(f"Total number of jobs executed: {total_executed}")
    if print_bool:
        print(f"Total number of jobs described: {total_printed}")


if __name__ == "__main__":
    main()
