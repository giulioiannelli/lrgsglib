from __future__ import annotations

from pathlib import Path

import subprocess

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    build_slanzarv_command,
    format_slanzarv_command,
)
from lrgsglib import *
from parsers.L3D_Recon import L3D_Recon_progName, L3D_Recon_progNameShrt
from parsers.L3D_Recon_srun import parse_arguments, parser


def _build_prog_args(L: int, pflip: float, temp: float) -> list[str]:
    return [str(L), f"{pflip:.3g}", f"{temp:.3g}"]


def main():
    args, unknown = parser.parse_known_args()
    exec_bool, print_bool = args.exec, args.print
    if not (exec_bool or print_bool):
        return

    progn = L3D_Recon_progName
    progn_shrt = L3D_Recon_progNameShrt

    side_list = list(args.dim_list)
    pflp_list = list(args.pflip_linsp)
    temp_list = list(args.Temp_linsp)

    memoryfunc = build_memory_function(args.slanzarv_minMB, args.slanzarv_maxMB, side_list)
    exec_path = LRGSG_SRC.relative_to(Path.cwd()) / f"{progn}.py"

    total_printed = total_executed = 0

    def dispatch(L: int, p: float, T: float) -> None:
        nonlocal total_printed, total_executed

        prog_args = _build_prog_args(L, p, T)
        cmd = ["python", str(exec_path), *prog_args, *unknown]

        slanz_opts = ["-m", str(memoryfunc(L))]
        if args.nomail:
            slanz_opts.append("--nomail")
        if args.short:
            slanz_opts.append("--short")
        if args.moretime:
            slanz_opts.extend(["--time", str(args.moretime)])
        jobname = build_jobname(
            program_short=progn_shrt,
            tokens=prog_args,
            job_id=args.slanzarv_id or None,
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


if __name__ == '__main__':
    main()
