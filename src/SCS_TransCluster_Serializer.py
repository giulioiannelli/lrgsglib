from __future__ import annotations

import shlex
import subprocess
from itertools import product
from pathlib import Path

import numpy as np

from parsers.SCS_TransCluster_Serializer import parser
from lrgsglib.config.progargs import SCS_TransCluster_progName, SCS_TransCluster_progNameShrt


def _memory_function(min_mb: int, max_mb: int, values: list[int]):
    if min_mb == max_mb:
        return lambda *_: min_mb

    bounds = [min(values), max(values)]
    targets = [min_mb, max_mb]

    def _interp(value: int) -> int:
        return int(np.interp(value, bounds, targets))

    return _interp


def _build_jobname(*, mode: str, N: int, gamma: float, j0: float, navg: int, suffix: str | None, job_id: str | None) -> str:
    tokens = [SCS_TransCluster_progNameShrt, mode, f"N{N}", f"g{gamma:.3g}", f"J0{j0:.3g}", f"navg{navg}"]
    if suffix:
        tokens.append(suffix)
    jobname = "_".join(tokens)
    if job_id:
        return f"{job_id}_{jobname}"
    return jobname


def main() -> None:
    args, unknown = parser.parse_known_args()

    if not (args.exec or args.print):
        return

    N_list = list(args.N_list)
    gamma_list = list(args.gamma_list)
    j0_list = list(args.J0_list)

    memoryfunc = _memory_function(args.slanzarv_minMB, args.slanzarv_maxMB, N_list)

    printed = executed = 0

    for N, gamma, j0 in product(N_list, gamma_list, j0_list):
        base_cmd = [
            "python",
            str(Path("src") / f"{SCS_TransCluster_progName}.py"),
            str(N),
            "--gamma",
            f"{gamma:.12g}",
            "--J0",
            f"{j0:.12g}",
            "--J",
            f"{args.J:.12g}",
            "--g",
            f"{args.g:.12g}",
            "--diagonal",
            args.diagonal,
            "-wd",
            args.workdir,
            "-na",
            str(args.number_of_averages),
            "-m",
            args.mode,
            "-sf",
            str(args.save_frequency),
            "-t",
            args.float_type,
            "--backend",
            args.backend,
            "--partition-rule",
            args.partition_rule,
        ]
        if args.out_suffix:
            base_cmd.extend(["-o", args.out_suffix])
        if args.seed is not None:
            base_cmd.extend(["--seed", str(args.seed)])
        base_cmd.extend(unknown)

        jobname = _build_jobname(
            mode=args.mode,
            N=N,
            gamma=gamma,
            j0=j0,
            navg=args.number_of_averages,
            suffix=args.out_suffix or None,
            job_id=args.slanzarv_id or None,
        )

        slanz_cmd = [
            "slanzarv",
            "-m",
            str(memoryfunc(N)),
        ]
        if args.nomail:
            slanz_cmd.append("--nomail")
        if args.short:
            slanz_cmd.append("--short")
        if args.moretime:
            slanz_cmd.extend(["--time", str(args.moretime)])
        slanz_cmd.extend(["--jobname", jobname])
        slanz_cmd.extend(base_cmd)

        if args.print:
            print(shlex.join(slanz_cmd))
            printed += 1
        if args.exec:
            subprocess.run(slanz_cmd, check=True)
            executed += 1

    if args.print:
        print(f"Total jobs described: {printed}")
    if args.exec:
        print(f"Total jobs executed: {executed}")


if __name__ == "__main__":
    main()
