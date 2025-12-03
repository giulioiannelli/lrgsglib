"""Kernel utilities for the generalized SCS transient cluster simulations."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

from lrgsglib import SCSGeneralizedNN
from lrgsglib.utils.basic.paths import remove_if_exists
from lrgsglib.utils.tools import ConditionalPartitioning

from .SCS import (
    SCSParameters,
    build_scs_kwargs,
    probe_output_graph,
)
from .generic import resolve_backend, resolve_float_type, find_existing_averages

__all__ = ["run_transcluster"]


def _compute_order_parameters(
    params: SCSParameters,
    *,
    backend: str,
    float_type: type,
    partitioner: ConditionalPartitioning,
    seed: int | None = None,
) -> tuple[float, float, float]:
    # build_scs_kwargs will use the provided seed if given, otherwise params.seed
    kwargs = build_scs_kwargs(params, make_dir_tree=False, seed=seed)
    scs = SCSGeneralizedNN(**kwargs)
    scs.compute_laplacian_spectrum_weigV(backend=backend, typf=float_type)
    scs.load_eigV_on_graph(which=0, binarize=False)
    scs.make_eigVclustersYN(val=partitioner, which=0, binarize=True)
    scs.compute_pinf(which=0, val=partitioner)

    gap = scs.get_gap()
    cluster_fraction = scs.Pinf
    energy0 = scs.get_sksph_energy_eigV(0)
    energy1 = scs.get_sksph_energy_eigV(1)
    return gap, cluster_fraction, energy0 - energy1


def _write_results(path: Path, data: Iterable[Iterable[float]]) -> None:
    array = np.asarray(list(data), dtype=float)
    if array.size == 0:
        return
    np.savetxt(path, array, fmt="%.7g")


def _make_build_output_path(
    base_path: Path, suffix_tokens: list[str], out_suffix: str | None
) -> Callable[[int], Path]:
    """Return a closure that builds output Path objects for a given na."""
    def build_output_path(na: int) -> Path:
        tokens = list(suffix_tokens) + [f"na={na}"]
        if out_suffix:
            tokens.append(out_suffix)
        return base_path / ("_".join(tokens) + ".txt")

    return build_output_path


def _read_prev_partial(prev_path: Path) -> tuple[np.ndarray, int]:
    """Load a previous partial file (new averages-only format) and return
    (accumulator_avg, na).

    New file layout (per refactor) is:
      [avg_smax, avg_smax2, std, avg_gap, avg_ediff]

    We extract na from the filename (expects a token 'na=<N>' in the stem)
    because na is no longer saved inside the file.
    """
    arr = np.loadtxt(prev_path).reshape(-1)
    # arr layout: [avg_smax, avg_smax2, std, avg_gap, avg_ediff]
    if arr.size < 5:
        raise ValueError(f"Previous partial file {prev_path} has unexpected format")

    # Build accumulator in our internal order: [smax_avg, smax2_avg, gap_avg, ediff_avg]
    accumulator_avg = np.array([
        float(arr[0]),  # avg_smax
        float(arr[1]),  # avg_smax2
        float(arr[3]),  # avg_gap
        float(arr[4]),  # avg_ediff
    ], dtype=float)

    # parse na from filename stem
    stem = prev_path.stem
    na = 0
    for part in stem.split('_'):
        if part.startswith('na='):
            try:
                na = int(part.split('=', 1)[1])
            except Exception:
                pass
            break

    if na <= 0:
        raise ValueError(f"Could not determine na from filename: {prev_path.name}")

    return accumulator_avg, na


def _compute_averages_and_stats(
    accumulator_avg: np.ndarray,
) -> tuple[float, float, float, float, float]:
    """Compute averages and std from normalized accumulator array.

    accumulator_avg is [smax_avg, smax2_avg, gap_avg, ediff_avg].
    Returns: (smax_avg, smax2_avg, std, gap_avg, ediff_avg)
    """
    smax_avg = float(accumulator_avg[0])
    smax2_avg = float(accumulator_avg[1])
    gap_avg = float(accumulator_avg[2])
    ediff_avg = float(accumulator_avg[3])
    var = smax2_avg - smax_avg**2
    std = np.sqrt(max(var, 0.0))
    return smax_avg, smax2_avg, std, gap_avg, ediff_avg


def _pack_results(
    smax_avg: float,
    smax2_avg: float,
    std: float,
    gap_avg: float,
    ediff_avg: float,
) -> list[list[float]]:
        """Return results in the nested-list format used by _write_results.

        New layout (averages only):
            [avg_smax, avg_smax2, std, avg_gap, avg_ediff]
        """
        return [[smax_avg, smax2_avg, std, gap_avg, ediff_avg]]


def _compute_and_write(path: Path, accumulator_avg: np.ndarray) -> None:
        """Compute averages from normalized accumulator array, pack and write to path.

        The written file contains averages only in the order:
            [avg_smax, avg_smax2, std, avg_gap, avg_ediff]
        """
        (
            smax_avg,
            smax2_avg,
            std,
            gap_avg,
            ediff_avg,
        ) = _compute_averages_and_stats(accumulator_avg)
        rows = _pack_results(smax_avg, smax2_avg, std, gap_avg, ediff_avg)
        _write_results(path, rows)



def _gen_sample_seed(rng: np.random.Generator | None) -> int | None:
    """Generate a random seed for sampling, or None if no RNG is provided."""
    if rng is None:
        return None
    return int(rng.integers(0, 2**32))


def _setup_output_paths(
    args, params: SCSParameters
) -> tuple[Callable[[int], Path], int]:
    """Setup output directory structure and determine existing progress.

    Returns:
        tuple: (build_output_path_func, existing_na)
    """
    probe = probe_output_graph(params)
    base_dir = probe.path_phtra
    if getattr(args, "verbose", False):
        print(f"Output directory: {base_dir}")

    suffix_tokens = [
        args.mode,
        probe.std_fname,
        f"J={args.J:.3g}",
        f"g={args.g:.3g}",
        f"J0={args.J0:.3g}",
    ]
    del probe
    base_filename = "_".join(suffix_tokens)

    base_path = Path(base_dir)

    # Build pattern for finding existing files
    pattern = f"{base_filename}_na=*"
    if args.out_suffix:
        pattern += f"_{args.out_suffix}"
    pattern += ".txt"

    existing_na = find_existing_averages(
        base_path,
        pattern,
        count_token="na",
        find_max=True,
        logger=None,
    )

    build_output_path = _make_build_output_path(
        base_path, suffix_tokens, args.out_suffix
    )

    return build_output_path, existing_na


def _initialize_accumulator(
    start_avg: int, build_output_path: Callable[[int], Path]
) -> np.ndarray:
    """Initialize accumulator from previous partial file if it exists.

    Returns:
        np.ndarray: accumulator array [smax_avg, smax2_avg, gap_avg, ediff_avg]
    """
    accumulator = np.zeros(4, dtype=float)
    if start_avg > 0:
        prev_path = build_output_path(start_avg)
        if prev_path.exists():
            accumulator, _ = _read_prev_partial(prev_path)
    return accumulator


def _update_running_averages(
    accumulator: np.ndarray,
    n: int,
    smax: float,
    gap: float,
    energy_diff: float,
) -> None:
    """Update running averages using incremental formula (in-place).

    accumulator layout: [smax_avg, smax2_avg, gap_avg, ediff_avg]
    """
    accumulator[0] += (smax - accumulator[0]) / n
    accumulator[1] += (smax * smax - accumulator[1]) / n
    accumulator[2] += (gap - accumulator[2]) / n
    accumulator[3] += (energy_diff - accumulator[3]) / n


def _save_intermediate_if_needed(
    args,
    avg_idx: int,
    current_total: int,
    accumulator: np.ndarray,
    build_output_path: Callable[[int], Path],
) -> None:
    """Save intermediate results and clean up old files if save frequency is met."""
    if not args.save_frequency or args.save_frequency <= 0:
        return

    if (avg_idx + 1) % args.save_frequency == 0:
        interim_path = build_output_path(current_total)
        _compute_and_write(interim_path, accumulator)

        # Remove older partial file
        old_n = current_total - args.save_frequency
        if old_n > 0:
            remove_if_exists(build_output_path(old_n))


def run_transcluster(args) -> None:
    """Main entry point for transient cluster order parameter calculation.

    Computes averaged order parameters (cluster size, gap, energy
    difference) for SCS generalized networks with incremental saving
    and resume capability.
    """
    # Validate inputs
    if args.mode != "ordParam":
        raise ValueError(
            "SCS_TransCluster currently supports only 'ordParam' mode."
        )
    if args.number_of_averages <= 0:
        raise ValueError("Number of averages must be positive.")

    # Setup computation parameters
    backend = resolve_backend(args.backend)
    float_type = resolve_float_type(args.float_type)
    partitioner = ConditionalPartitioning(args.partition_rule)

    params = SCSParameters(
        N=args.N,
        J0=args.J0,
        gamma=args.gamma,
        J=args.J,
        g=args.g,
        diagonal=args.diagonal,
        workdir=args.workdir,
        seed=args.seed,
    )

    # Setup output paths and check for existing results
    build_output_path, existing_na = _setup_output_paths(args, params)

    if existing_na >= args.number_of_averages:
        if getattr(args, "verbose", False):
            print(
                f"File with {existing_na} averages already exists. "
                "Nothing to do."
            )
        return

    start_avg = existing_na
    remaining_averages = args.number_of_averages - existing_na
    if existing_na > 0 and getattr(args, "verbose", False):
        print(
            f"Found existing file with {existing_na} averages. "
            f"Computing {remaining_averages} more..."
        )

    final_output_path = build_output_path(args.number_of_averages)
    if final_output_path.exists() and not getattr(
        args, "remove_files", False
    ):
        raise SystemExit(f"File {final_output_path.name} already exists.")

    # Initialize accumulator (resume from previous if available)
    accumulator = _initialize_accumulator(start_avg, build_output_path)

    # Main computation loop
    rng = (
        np.random.default_rng(args.seed)
        if args.seed is not None
        else None
    )

    for avg_idx in range(remaining_averages):
        sample_seed = _gen_sample_seed(rng)

        gap, smax, energy_diff = _compute_order_parameters(
            params,
            backend=backend,
            float_type=float_type,
            partitioner=partitioner,
            seed=sample_seed,
        )

        # Update running averages incrementally
        current_total = start_avg + avg_idx + 1
        _update_running_averages(
            accumulator, current_total, smax, gap, energy_diff
        )

        # Save intermediate results periodically
        _save_intermediate_if_needed(
            args, avg_idx, current_total, accumulator, build_output_path
        )

    # Write final results
    _compute_and_write(final_output_path, accumulator)

    # Clean up old partial file
    if start_avg > 0:
        remove_if_exists(build_output_path(start_avg))

