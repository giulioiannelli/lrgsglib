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
    auto_select_gpu_type,
)
from parsers.MCG_SlaplSpect_Serializer import (
    MCG_MCGSS_progname,
    MCG_MCGSS_progname_shrt,
    parser,
)
# Import directly from config modules to avoid loading entire lrgsglib package
from lrgsglib.config.progargs.SlaplSpect import (
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
    # Filter out the '--' separator if present
    unknown = [arg for arg in unknown if arg != '--']

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

    script_path = Path("src") / f"{MCG_MCGSS_progname}.py"

    total_printed = total_executed = 0

    def estimate_memory_and_backend(p1: float, p2: float, p3: float, p4: float,
                                     fraction: float, iterations: int) -> tuple[int, str, bool]:
        """
        Estimate memory requirement and determine optimal backend.

        Uses analytical formulas with giant component correction:
        - SIZE = 2^(iterations + 1)
        - N_grid = fraction × SIZE²
        - N_actual ≈ N_grid for fraction ≥ 0.6 (giant component ≈ full graph)
        - N_actual << N_grid for fraction < 0.6 (but already small, so conservative overestimate is safe)

        Returns
        -------
        tuple[int, str, bool]
            (memory_mb, backend, use_gpu_node)
            - memory_mb: RAM to allocate
            - backend: 'cupy' or 'scipy'
            - use_gpu_node: whether to request GPU node from slanzarv
        """
        # Compute graph size analytically
        SIZE = 2 ** (iterations + 1)
        N_grid = max(1, round(fraction * SIZE * SIZE))

        # Apply giant component correction
        # For fraction ≥ 0.6: giant component ≈ 95-100% of selected nodes
        # For fraction < 0.6: giant component < 60%, but already small so overestimate is safe
        #
        # P-matrix impact (measured at fraction=0.7, it=7):
        # - User's matrices (p1=1.0, p2=p3∈[0.7,0.9], p4=0.9): N ∈ [42.7k, 44.8k], ~3% variation
        # - Overall variation across realistic matrices: ~10% coefficient of variation
        # - Safety margin: 1.15× accounts for p-matrix uncertainty
        if fraction >= 0.6:
            # Apply 15% safety margin for p-matrix variability
            N = int(N_grid * 1.15)
        else:
            # Conservative overestimate (actual will be much smaller, but allocate for full size)
            N = int(N_grid * 1.15)

        E_approx = 2 * N  # Approximate for grid-based graphs

        # Calculate matrix sizes
        dense_gb = (N ** 2 * 8) / (1024 ** 3)
        sparse_mb = (E_approx * 16) / (1024 ** 2)

        # Apply same thresholds as select_optimal_backend
        # A100-80GB: Use 70GB limit (leave 10GB for workspace/overhead)
        # Dense matrix + workspace ≈ 3× matrix size, so N_max ≈ sqrt(70/3) ≈ 48k nodes
        GPU_VRAM_LIMIT_GB = 70.0
        # High-memory CPU nodes (Nix/Orion) have 768-1536 GB RAM
        CPU_RAM_LIMIT_GB = 500.0

        # NOTE: This is for FULL SPECTRUM (--howmany 0)
        # For full spectrum, dense is ALWAYS better than sparse
        # (sparse with k≈N allocates N×N workspace)

        if dense_gb < GPU_VRAM_LIMIT_GB:
            # Dense GPU: CPU builds graph, GPU computes eigenvalues
            # Allocate 2x dense matrix size for graph construction + workspace
            memory_mb = int(dense_gb * 1024 * 2.0)
            return (max(memory_mb, 2048), 'cupy', True)  # Request GPU node
        elif dense_gb < CPU_RAM_LIMIT_GB:
            # Dense CPU: matrix + 3.5x workspace for eigensolver
            # High-memory nodes (Orion: 1536GB) can handle this
            # For full spectrum, dense is always better than sparse!
            memory_mb = int(dense_gb * 1024 * 3.5)  # 3.5x for safety
            return (min(memory_mb, 1200 * 1024), 'scipy', False)  # CPU-only node
        else:
            # Graph too large for full spectrum
            # User should use partial spectrum (--howmany) instead
            memory_mb = 300 * 1024  # Placeholder, will fail
            return (memory_mb, 'scipy', False)

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

        # Check if user specified backend in unknown args
        user_backend = None
        for i, arg in enumerate(unknown):
            if arg == '--backend' and i + 1 < len(unknown):
                user_backend = unknown[i + 1]
                break

        # Determine optimal backend and memory allocation
        memory_mb, backend, use_gpu_node = estimate_memory_and_backend(
            p1, p2, p3, p4, fraction, iterations
        )

        # If user specified backend, override the determined backend
        if user_backend is not None:
            backend = user_backend
            # Recalculate use_gpu_node based on user's backend choice
            use_gpu_node = (backend == 'cupy')

            # Recalculate memory with user's backend choice
            memory_mb, _, _ = estimate_memory_and_backend(
                p1, p2, p3, p4, fraction, iterations
            )

        # Build explicit command with all options visible
        # Start with backend (determined or user-specified)
        explicit_opts = ["--backend", backend]

        # Parse unknown args to extract user-specified options
        user_opts = {}
        i = 0
        while i < len(unknown):
            if unknown[i].startswith('--'):
                key = unknown[i]
                # Check if next arg is a value (not another flag)
                if i + 1 < len(unknown) and not unknown[i + 1].startswith('--'):
                    user_opts[key] = unknown[i + 1]
                    i += 2
                else:
                    # Boolean flag or flag without value
                    user_opts[key] = True
                    i += 1
            else:
                i += 1

        # Add user-specified options explicitly (skip --backend if already added)
        for key, value in user_opts.items():
            if key != '--backend':
                if value is True:
                    explicit_opts.append(key)
                else:
                    explicit_opts.extend([key, str(value)])

        cmd = ["python", str(script_path), *prog_args, *explicit_opts]

        final_cmd = cmd
        if use_slanzarv:
            slanz_opts = ["-m", str(memory_mb)]
            if args.nomail:
                slanz_opts.append("--nomail")
            if args.short:
                slanz_opts.append("--short")
            if args.moretime:
                slanz_opts.extend(["--time", str(args.moretime)])

            # Determine partition: AMD for high-memory jobs (>100GB)
            # User can override with --partition flag
            partition = args.partition
            if partition is None and memory_mb > 100 * 1024:  # 100GB in MB
                partition = 'AMD'

            if partition is not None:
                slanz_opts.extend(["--partition", partition])

            # Automatically request GPU node when backend='cupy'
            # Backend is determined by graph size or user override via --backend
            if use_gpu_node and backend == 'cupy':
                slanz_opts.append("--gpu")

                # Select GPU type
                if args.gputype:
                    gputype = args.gputype
                else:
                    # Default to A100 (compatible with CUDA 12.x in current conda environment)
                    gputype = 'a100'

                slanz_opts.extend(["--gputype", gputype])
            else:
                # CPU-only job: specify core count based on job size
                # Smaller jobs (< 50GB): 4 cores
                # Larger jobs (>= 50GB): 8 cores
                cores = 4 if memory_mb < 50 * 1024 else 8
                slanz_opts.extend(["-c", str(cores)])

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
