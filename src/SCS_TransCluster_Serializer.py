from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

from kernels.Serializer import build_jobname, build_memory_function
from parsers.SCS_TransCluster_Serializer import parser, SCS_TransCluster_progName, SCS_TransCluster_progNameShrt
from lrgsglib.config.progargs import (
    DEFAULT_SCS_NN_SRUN_N_LIST,
    DEFAULT_SCS_NN_SRUN_GAMMA_LIST,
    DEFAULT_SCS_NN_SRUN_J0_LIST,
    DEFAULT_SCS_NN_SRUN_J_LIST,
    DEFAULT_SCS_NN_SRUN_G_LIST,
    DEFAULT_SCS_NN_SRUN_DIAGONAL_LIST,
)


def _collect_values(
    explicit: Iterable[float] | None,
    linspace_values,
    default_values: Sequence[float],
) -> list[float]:
    return _collect_values_typed(explicit, linspace_values, default_values, float)


def _collect_values_typed(
    explicit: Iterable[float] | None,
    linspace_values,
    default_values: Sequence[float],
    caster: Callable[[float], float],
) -> list[float]:
    values: list[float] = []
    if explicit:
        values.extend(caster(v) for v in explicit)
    if linspace_values is not None:
        values.extend(caster(v) for v in np.asarray(linspace_values, dtype=float))
    if not values:
        values.extend(caster(v) for v in default_values)
    return values


def main() -> None:
    args, unknown = parser.parse_known_args()
    exec_bool, print_bool = args.exec, args.print
    if not (exec_bool or print_bool):
        return

    N_values = _collect_values_typed(
        args.N_list,
        None,  # No linspace for N
        DEFAULT_SCS_NN_SRUN_N_LIST,
        int,
    )

    gamma_values = _collect_values_typed(
        args.gamma_list,
        None,  # No linspace for gamma
        DEFAULT_SCS_NN_SRUN_GAMMA_LIST,
        float,
    )

    j0_values = _collect_values(
        args.J0_list,
        None,  # No linspace for J0
        DEFAULT_SCS_NN_SRUN_J0_LIST,
    )

    J_values = _collect_values(
        args.J_list,
        args.J_linsp,
        DEFAULT_SCS_NN_SRUN_J_LIST,
    )

    g_values = _collect_values(
        args.g_list,
        args.g_linsp,
        DEFAULT_SCS_NN_SRUN_G_LIST,
    )

    # Diagonal values are strings
    diagonal_values: list[str] = []
    if args.diagonal_list:
        diagonal_values.extend(str(v) for v in args.diagonal_list)
    if not diagonal_values:
        diagonal_values.extend(str(v) for v in DEFAULT_SCS_NN_SRUN_DIAGONAL_LIST)

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
            f"{gamma:.12g}",
            "--J0",
            f"{j0:.12g}",
            "--J",
            f"{J_value:.12g}",
            "--g",
            f"{g_value:.12g}",
            "--diagonal",
            diagonal,
        ]

        cmd = ["python", str(script_path), *prog_args, *unknown]

        slanz_opts = ["-m", str(memoryfunc(N))]
        if args.nomail:
            slanz_opts.append("--nomail")
        if args.short:
            slanz_opts.append("--short")
        if args.moretime:
            slanz_opts.extend(["--time", str(args.moretime)])

        # Build jobname from N, gamma, J, g, j0, diagonal
        jobname_tokens = [
            f"N{N}",
            f"g{gamma:.3g}",
            f"J{J_value:.3g}",
            f"gcpl{g_value:.3g}",
            f"J0{j0:.3g}",
            f"diag{diagonal}",
        ]
        jobname = build_jobname(
            program_short=SCS_TransCluster_progNameShrt,
            tokens=jobname_tokens,
            job_id=args.slanzarv_id or None,
            prefix_job_id=True,
        )
        slanz_opts.extend(["--jobname", jobname])

        slanz_cmd = ["slanzarv", *slanz_opts, *cmd]

        if print_bool:
            print(" ".join(slanz_cmd))
            total_printed += 1
        if exec_bool:
            subprocess.run(slanz_cmd, check=True)
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
