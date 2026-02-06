"""Shared helpers for serializer-style launcher scripts."""

from __future__ import annotations

from typing import Callable, Iterable, Sequence

import numpy as np

__all__ = [
    "build_memory_function",
    "build_jobname",
    "resolve_slanzarv_mode",
    "_determine_precision",
    "_format_value_consistently",
    "_collect_values",
    "_collect_values_typed",
    "format_slanzarv_command",
    "build_slanzarv_command",
    "estimate_mcg_laplacian_memory_mb",
    "auto_select_gpu_type",
]

def build_memory_function(min_mb: int, max_mb: int, values: Sequence[int]) -> Callable[[int], int]:
    """
    Return an interpolating function that maps parameter values to memory limits.

    When the provided min/max limits are identical, a constant function is returned.
    """
    if min_mb == max_mb or not values:
        return lambda *_: min_mb

    low, high = min(values), max(values)
    if low == high:
        return lambda *_: min_mb

    bounds = [low, high]
    targets = [min_mb, max_mb]

    def _interp(value: int) -> int:
        return int(np.interp(value, bounds, targets))

    return _interp


def build_jobname(
    *,
    program_short: str,
    tokens: Iterable[str],
    job_id: str | None = None,
    delimiter: str = "_",
    prefix_job_id: bool = False,
) -> str:
    """
    Compose a Slurm job name using a program short name and additional tokens.

    Set ``prefix_job_id`` to True if the external job identifier should appear
    before the program short name (matching historical behaviour in some scripts).
    """
    cleaned_tokens = [tok for tok in tokens if tok]
    parts = [program_short, *cleaned_tokens]

    if job_id:
        if prefix_job_id:
            parts = [job_id, *parts]
        else:
            parts = [program_short, job_id, *cleaned_tokens]

    return delimiter.join(parts)


def resolve_slanzarv_mode(
    mode: str | None,
    *,
    prefix: str = "slanzarv_",
    default_use_slanzarv: bool = False,
) -> tuple[bool, str | None]:
    """
    Determine whether slanzarv should be used based on the ``mode`` argument.

    Parameters
    ----------
    mode:
        The user-specified mode string (may include the ``slanzarv_`` prefix).
    prefix:
        Prefix that triggers slanzarv submission.
    default_use_slanzarv:
        Whether to default to slanzarv when ``mode`` is missing.

    Returns
    -------
    tuple[bool, str | None]
        ``(use_slanzarv, stripped_mode)`` where ``stripped_mode`` has the prefix
        removed (or is ``None`` if unavailable).
    """
    if not mode:
        return default_use_slanzarv, None

    if mode.startswith(prefix):
        stripped = mode[len(prefix):] or None
        return True, stripped

    return False, mode


def _determine_precision(values: list[float]) -> int:
    """
    Determine the minimum number of decimal places needed to distinguish all values.
    Returns the number of decimal places needed.
    """
    if not values:
        return 2

    # Determine intrinsic precision from the values themselves
    # by checking how many decimal places are present in the string representation
    intrinsic_precision = 2
    for val in values:
        # Convert to string and count decimal places
        val_str = f"{val:.15g}"  # Use general format to avoid trailing zeros
        if '.' in val_str:
            decimal_part = val_str.split('.')[1]
            # Remove scientific notation if present
            if 'e' in decimal_part.lower():
                decimal_part = decimal_part.split('e')[0]
            intrinsic_precision = max(intrinsic_precision, len(decimal_part))

    if len(values) <= 1:
        # For single value, use its intrinsic precision
        return min(intrinsic_precision, 8)

    # Sort values to find minimum differences
    sorted_vals = sorted(set(values))
    if len(sorted_vals) == 1:
        return min(intrinsic_precision, 8)

    # Find the minimum difference between consecutive values
    min_diff = min(abs(sorted_vals[i+1] - sorted_vals[i]) for i in range(len(sorted_vals) - 1))

    # Determine precision needed: we need enough decimals so that min_diff is distinguishable
    if min_diff == 0:
        return min(intrinsic_precision, 8)

    # Calculate how many decimal places we need
    # We want at least 2 significant digits in the minimum difference
    precision = max(2, int(np.ceil(-np.log10(min_diff) + 1)))

    # Use the maximum of intrinsic precision and distinction precision
    precision = max(precision, intrinsic_precision)

    # Cap at reasonable maximum
    return min(precision, 8)


def _format_value_consistently(value: float, precision: int) -> str:
    """
    Format a value with consistent decimal places.
    All values use exactly the same number of decimal places (the precision).
    """
    return f"{value:.{precision}f}"


def _collect_values(
    explicit: Iterable[float] | None,
    linspace_values,
    default_values: Sequence[float] | None = None,
) -> list[float]:
    """Collect numeric values from explicit inputs or a linspace.

    Falls back to ``default_values`` when neither explicit nor linspace inputs
    are provided. Mirrors :func:`_collect_values_typed` but defaults to
    ``float`` casting.
    """

    return _collect_values_typed(explicit, linspace_values, default_values or [], float)


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


def build_slanzarv_command(
    slanz_opts: list[str],
    cmd: list[str],
) -> list[str]:
    """
    Build a complete slanzarv command with sbatch options.

    Parameters
    ----------
    slanz_opts : list[str]
        Slanzarv-specific options (e.g., ["-m", "2048", "--nomail", "--jobname", "myjob"])
    cmd : list[str]
        The actual command to execute (e.g., ["python", "script.py", "arg1", "arg2"])

    Returns
    -------
    list[str]
        Complete command list including slanzarv, sbatch options, and the command
    """
    sbatch_opts = [
        "--output=.log/%x_%j.out",
        "--error=.log/%x_%j.err",
    ]

    return ["slanzarv", *slanz_opts, *sbatch_opts, "--", *cmd]


def format_slanzarv_command(
    slanz_opts: list[str],
    cmd: list[str],
) -> str:
    """
    Format a slanzarv command for readable multi-line output.

    Produces output in the format:
        slanzarv <options> --output=... --error=... -- \\
          python script.py args...

    Parameters
    ----------
    slanz_opts : list[str]
        Slanzarv-specific options (e.g., ["-m", "2048", "--nomail", "--jobname", "myjob"])
    cmd : list[str]
        The actual command to execute (e.g., ["python", "script.py", "arg1", "arg2"])

    Returns
    -------
    str
        Formatted multi-line command string with backslash continuation
    """
    sbatch_opts = [
        "--output=.log/%x_%j.out",
        "--error=.log/%x_%j.err",
    ]

    slanzarv_sbatch_line = " ".join(["slanzarv", *slanz_opts, *sbatch_opts, "--"])
    python_line = "  " + " ".join(cmd)

    return f"{slanzarv_sbatch_line} \\\n{python_line}"


def estimate_mcg_laplacian_memory_mb(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    iterations: int,
    fraction: float = 0.4,
    stochastic: bool = False,
) -> float:
    """
    Estimate memory requirements (in MB) for MCG Laplacian eigenvalue computation.

    Generates a single sample graph to get actual N, then estimates memory based on
    dense Laplacian matrix storage + eigenvalue solver workspace.

    Parameters
    ----------
    p1, p2, p3, p4 : float
        Cascade probability matrix parameters
    iterations : int
        Number of cascade iterations
    fraction : float
        Fraction of edges kept (default 0.4)
    stochastic : bool
        Whether to use stochastic mode (default False for deterministic sizing)

    Returns
    -------
    float
        Estimated memory in MB (conservative upper bound)

    Notes
    -----
    Memory includes:
    - Dense Laplacian matrix: N × N × 8 bytes (float64)
    - Eigenvalue solver workspace: ~3-4× matrix size (for full decomposition)
    - Safety margin: 20%

    This function actually generates a sample graph, so it's accurate but may take
    a few seconds for large iterations (>10).
    """
    from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

    # Generate a sample graph to get actual N
    mc = MultiplicativeCascadeGraph(
        p1=p1,
        p2=p2,
        p3=p3,
        p4=p4,
        fraction=fraction,
        iterations=iterations,
        stochastic=stochastic,
    )

    N = mc.N

    # Memory for dense Laplacian matrix (N x N, float64 = 8 bytes)
    matrix_mb = (N ** 2 * 8) / (1024 ** 2)

    # Eigenvalue solver workspace (conservative: 4x matrix size for full decomposition)
    # This accounts for intermediate arrays, eigenvectors, working memory, etc.
    total_mb = matrix_mb * 4

    # Add safety margin (20%)
    total_mb *= 1.2

    return total_mb


def auto_select_gpu_type(
    p1: float,
    p2: float,
    p3: float,
    p4: float,
    iterations: int,
    fraction: float = 0.4,
    stochastic: bool = False,
    fermi_threshold_mb: float = 1024.0,
    verbose: bool = False
) -> str:
    """
    Automatically select GPU type based on estimated memory requirements.

    Generates a sample graph to get actual size, then selects appropriate GPU.

    Parameters
    ----------
    p1, p2, p3, p4 : float
        Cascade probability matrix parameters
    iterations : int
        Number of cascade iterations
    fraction : float
        Fraction of edges kept (default 0.4)
    stochastic : bool
        Whether to use stochastic mode (default False for deterministic sizing)
    fermi_threshold_mb : float
        Memory threshold (MB) for selecting Fermi vs A100 (default 1024 = 1GB)
    verbose : bool
        Print selection reasoning

    Returns
    -------
    str
        GPU type: 'fermi' for small problems, 'a100' for large problems

    Notes
    -----
    Fermi (Tesla C2050):
    - 3GB VRAM, CUDA 8.0 max
    - Suitable for small/medium graphs (< 1GB memory)

    A100:
    - 80GB VRAM, CUDA 12.1+
    - Needed for large graphs or when using modern CUDA/CuPy

    Conservative threshold of 1GB ensures Fermi isn't overloaded.
    """
    estimated_mb = estimate_mcg_laplacian_memory_mb(
        p1=p1,
        p2=p2,
        p3=p3,
        p4=p4,
        iterations=iterations,
        fraction=fraction,
        stochastic=stochastic,
    )

    if estimated_mb < fermi_threshold_mb:
        gpu_type = 'fermi'
        reason = f"small problem (estimated {estimated_mb:.1f}MB < {fermi_threshold_mb:.0f}MB threshold)"
    else:
        gpu_type = 'a100'
        reason = f"large problem (estimated {estimated_mb:.1f}MB >= {fermi_threshold_mb:.0f}MB threshold)"

    if verbose:
        print(f"Auto-selecting GPU: {gpu_type} ({reason})")

    return gpu_type
