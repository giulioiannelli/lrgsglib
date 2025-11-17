"""Kernel utilities for the generalized SCS transient cluster simulations."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

from lrgsglib import SCSGeneralizedNN
from lrgsglib.utils.tools import ConditionalPartitioning

from .SCS import (
    SCSParameters,
    build_scs_kwargs,
    probe_output_graph,
)
from .generic import resolve_backend, resolve_float_type

__all__ = ["run_transcluster"]


def _largest_cluster_fraction(graph: SCSGeneralizedNN) -> float:
    clusters: Sequence[Sequence[int]] = getattr(graph, "clustersY", [])
    if not clusters:
        return 0.0
    largest = max(clusters, key=len)
    return len(largest) / float(graph.N)


def _compute_order_parameters(
    params: SCSParameters,
    *,
    j0: float,
    backend: str,
    float_type: type,
    partitioner: ConditionalPartitioning,
    seed: int | None,
) -> tuple[float, float, float]:
    kwargs = build_scs_kwargs(params, J0_value=j0, make_dir_tree=False, seed=seed)
    scs = SCSGeneralizedNN(**kwargs)
    scs.compute_laplacian_spectrum_weigV(backend=backend, typf=float_type)
    scs.load_eigV_on_graph(which=0, binarize=False)
    scs.make_eigVclustersYN(val=partitioner, which=0, binarize=True)

    gap = (scs.eigv[1] - scs.eigv[0]) / scs.eigv[-1]
    cluster_fraction = _largest_cluster_fraction(scs)
    energy0 = scs.get_sksph_energy_eigV(0)
    energy1 = scs.get_sksph_energy_eigV(1)
    return gap, cluster_fraction, energy0 - energy1


def _write_results(path: Path, data: Iterable[Iterable[float]]) -> None:
    array = np.asarray(list(data), dtype=float)
    if array.size == 0:
        return
    np.savetxt(path, array, fmt="%.7g")


def _find_existing_navg(base_path: Path, base_filename: str, out_suffix: str | None) -> int:
    """Search for existing partial result files and return the largest found navg value."""
    existing_navg = 0
    if not base_path.exists():
        return 0

    import glob

    pattern = str(base_path / f"{base_filename}_navg=*")
    if out_suffix:
        pattern += f"_{out_suffix}"
    pattern += ".txt"

    existing_files = glob.glob(pattern)
    for filepath in existing_files:
        filename = Path(filepath).stem
        try:
            navg_part = [p for p in filename.split('_') if p.startswith("navg=")][0]
            navg_val = int(navg_part.split('=')[1])
            existing_navg = max(existing_navg, navg_val)
        except Exception:
            # ignore unparsable files
            pass

    return existing_navg


def _make_build_output_path(base_path: Path, suffix_tokens: list[str], out_suffix: str | None) -> Callable[[int], Path]:
    """Return a closure that builds output Path objects for a given navg."""
    def build_output_path(navg: int) -> Path:
        tokens = list(suffix_tokens) + [f"navg={navg}"]
        if out_suffix:
            tokens.append(out_suffix)
        return base_path / ("_".join(tokens) + ".txt")

    return build_output_path


def _read_prev_partial(prev_path: Path) -> tuple[np.ndarray, int]:
    """Load a previous partial file (new averages-only format) and return
    (accumulator_avg, navg).

    New file layout (per refactor) is:
      [avg_smax, avg_smax2, std, avg_gap, avg_ediff]

    We extract navg from the filename (expects a token 'navg=<N>' in the stem)
    because navg is no longer saved inside the file.
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

    # parse navg from filename stem
    stem = prev_path.stem
    navg = 0
    for part in stem.split('_'):
        if part.startswith('navg='):
            try:
                navg = int(part.split('=', 1)[1])
            except Exception:
                pass
            break

    if navg <= 0:
        raise ValueError(f"Could not determine navg from filename: {prev_path.name}")

    return accumulator_avg, navg


def _compute_averages_and_stats(accumulator_avg: np.ndarray) -> tuple[float, float, float, float, float]:
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


def _pack_results(smax_avg: float, smax2_avg: float, std: float, gap_avg: float, ediff_avg: float) -> list[list[float]]:
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
        smax_avg, smax2_avg, std, gap_avg, ediff_avg = _compute_averages_and_stats(
                accumulator_avg
        )
        rows = _pack_results(smax_avg, smax2_avg, std, gap_avg, ediff_avg)
        _write_results(path, rows)



def _gen_sample_seed(rng: np.random.Generator | None) -> int | None:
    if rng is None:
        return None
    return int(rng.integers(0, 2**32))


def _remove_if_exists(path: Path) -> None:
    if path.exists():
        path.unlink()


def run_transcluster(args) -> None:
    if args.mode != "ordParam":
        raise ValueError("SCS_TransCluster currently supports only 'ordParam' mode.")
    if args.number_of_averages <= 0:
        raise ValueError("Number of averages must be positive.")

    backend = resolve_backend(args.backend)
    float_type = resolve_float_type(args.float_type)
    partitioner = ConditionalPartitioning(args.partition_rule)
    j0_value = float(args.J0)

    params = SCSParameters(
        N=args.N,
        gamma=args.gamma,
        J=args.J,
        g=args.g,
        diagonal=args.diagonal,
        workdir=args.workdir,
        seed=args.seed,
    )

    probe = probe_output_graph(params, base_J0=j0_value)
    base_dir = probe.path_phtra

    suffix_tokens = [
        args.mode,
        probe.std_fname,
        f"J={args.J:.3g}",
        f"g={args.g:.3g}",
        f"J0={j0_value:.3g}",
    ]
    del probe
    base_filename = "_".join(suffix_tokens)

    # ---- READ EXISTING PARTIAL FILES ----
    base_path = Path(base_dir)

    existing_navg = _find_existing_navg(base_path, base_filename, args.out_suffix)

    if existing_navg >= args.number_of_averages:
        if getattr(args, "verbose", False):
            print(f"File with {existing_navg} averages already exists. Nothing to do.")
        return

    start_avg = existing_navg
    remaining_averages = args.number_of_averages - existing_navg
    if existing_navg > 0:
        if getattr(args, "verbose", False):
            print(f"Found existing file with {existing_navg} averages. Computing {remaining_averages} more...")

    build_output_path = _make_build_output_path(base_path, suffix_tokens, args.out_suffix)
    final_output_path = build_output_path(args.number_of_averages)
    if final_output_path.exists() and not getattr(args, "remove_files", False):
        raise SystemExit(f"File {final_output_path.name} already exists.")

    # ---------------------------------------------------------
    # INITIALIZE ACCUMULATORS, INCLUDING READING PREVIOUS DATA
    # ---------------------------------------------------------
    # accumulator now stores running averages: [smax_avg, smax2_avg, gap_avg, ediff_avg]
    accumulator = np.zeros(4, dtype=float)
    cluster_fractions = []

    # if previous partial file exists, read and initialize
    if start_avg > 0:
        prev_path = build_output_path(start_avg)
        if prev_path.exists():
            accumulator, prev_N = _read_prev_partial(prev_path)

    # ---------------------------------------------------------
    # MAIN LOOP
    # ---------------------------------------------------------
    rng = np.random.default_rng(args.seed) if args.seed is not None else None

    for avg_idx in range(remaining_averages):
        sample_seed = _gen_sample_seed(rng)

        gap, smax, energy_diff = _compute_order_parameters(
            params,
            j0=j0_value,
            backend=backend,
            float_type=float_type,
            partitioner=partitioner,
            seed=sample_seed,
        )

        # update running averages using incremental formula
        current_total_averages = start_avg + avg_idx + 1
        n = current_total_averages
        # accumulator holds previous average for n-1 samples (or 0). Update to include new sample
        # accumulator = [avg_smax, avg_smax2, avg_gap, avg_ediff]

        accumulator[0] += (smax - accumulator[0]) / n         # smax_avg
        accumulator[1] += (smax*smax - accumulator[1]) / n    # smax2_avg
        accumulator[2] += (gap - accumulator[2]) / n          # gap_avg
        accumulator[3] += (energy_diff - accumulator[3]) / n  # ediff_avg


        # ---- SAVING INTERMEDIATE ----
        if args.save_frequency and args.save_frequency > 0:
            if (avg_idx + 1) % args.save_frequency == 0:

                n_tot = current_total_averages
                interim_output_path = build_output_path(n_tot)
                _compute_and_write(interim_output_path, accumulator)

                # remove older partial file
                old_n = n_tot - args.save_frequency
                if old_n > 0:
                    old_path = build_output_path(old_n)
                    _remove_if_exists(old_path)

    # ---------------------------------------------------------
    # FINAL OUTPUT
    # ---------------------------------------------------------
    total_N = start_avg + remaining_averages

    _compute_and_write(final_output_path, accumulator)
    # remove old file with fewer averages
    if start_avg > 0:
        old_final_path = build_output_path(start_avg)
        _remove_if_exists(old_final_path)

