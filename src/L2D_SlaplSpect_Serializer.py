from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

from kernels.Serializer import (
    build_jobname,
    build_memory_function,
    resolve_slanzarv_mode,
    _collect_values,
    _collect_values_typed,
)
from parsers.L2D_SlaplSpect_Serializer import (
    parser,
    L2D_SlaplSpect_progName,
    L2D_SlaplSpect_progNameShrt,
)
from lrgsglib.config.progargs import (
    DEFAULT_L2D_SLAPLSPECT_SRUN_L_LIST,
    DEFAULT_L2D_SLAPLSPECT_SRUN_P_LIST,
)


# Per-mode default sweeps preserve the legacy serializer behaviour for
# `eigval_dist` and `eigvec_dist`. CLI lists (`--L-list`, `--L-linsp`,
# `--p-list`, `--p-linsp`) override these on a per-axis basis.
_DEFAULT_L_BY_MODE = {
    "eigvec_dist": [2 ** i for i in range(4, 10)],
    "eigval_dist": [16, 32, 48, 64, 96, 128],
}
_DEFAULT_P_BY_MODE = {
    "eigvec_dist": list(np.linspace(0.06, 0.115, 10)),
    "eigval_dist": [
        0.01, 0.025, 0.05, 0.075, 0.08, 0.09, 0.1, 0.11,
        0.125, 0.15, 0.25, 0.5, 0.75,
    ],
}


def main() -> None:
    args, unknown = parser.parse_known_args()
    if not (args.exec or args.print):
        return

    use_slanzarv, actual_mode = resolve_slanzarv_mode(
        getattr(args, "mode", None),
        default_use_slanzarv=False,
    )
    if not actual_mode:
        parser.error(
            "--mode is required (e.g. 'eigvals', 'eigval_dist', "
            "or prefixed with 'slanzarv_')"
        )

    default_L = _DEFAULT_L_BY_MODE.get(actual_mode, DEFAULT_L2D_SLAPLSPECT_SRUN_L_LIST)
    default_p = _DEFAULT_P_BY_MODE.get(actual_mode, DEFAULT_L2D_SLAPLSPECT_SRUN_P_LIST)

    L_values = _collect_values_typed(args.L_list, args.L_linsp, default_L, int)
    p_values = _collect_values(args.p_list, args.p_linsp, default_p)

    memoryfunc = build_memory_function(
        args.slanzarv_minMB, args.slanzarv_maxMB, L_values
    )
    script_path = Path("src") / f"{L2D_SlaplSpect_progName}.py"

    total_printed = total_executed = 0

    def dispatch(L: int, p: float) -> None:
        nonlocal total_printed, total_executed

        probability_str = f"{p:.3g}"
        prog_args = [
            str(L),
            "-p", probability_str,
            "--mode", actual_mode,
        ]
        cmd = ["python", str(script_path), *prog_args, *unknown]

        if use_slanzarv:
            slanz_opts = ["-m", str(memoryfunc(L))]
            if args.nomail:
                slanz_opts.append("--nomail")
            if args.short:
                slanz_opts.append("--short")
            if getattr(args, "moretime", 0):
                slanz_opts.extend(["--time", str(args.moretime)])

            jobname = build_jobname(
                program_short=L2D_SlaplSpect_progNameShrt,
                tokens=[actual_mode, f"L{L}", f"p{probability_str}"],
                job_id=getattr(args, "slanzarv_id", None) or None,
            )
            slanz_opts.extend(["--jobname", jobname])

            final_cmd = ["slanzarv", *slanz_opts, *cmd]
        else:
            final_cmd = cmd

        if args.print:
            print(" ".join(final_cmd))
            total_printed += 1
        if args.exec:
            subprocess.run(final_cmd, check=True)
            total_executed += 1

    for L in L_values:
        for p in p_values:
            dispatch(L, p)

    if args.exec:
        print(f"Total number of jobs executed: {total_executed}")
    if args.print:
        print(f"Total number of jobs described: {total_printed}")


if __name__ == "__main__":
    main()
