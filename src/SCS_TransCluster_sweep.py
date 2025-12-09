#!/usr/bin/env python3
"""CPU-friendly sweep runner for the generalized SCS transient cluster analysis.

This script mirrors the usability of ``test_cp2d_relu.py`` but targets the SCS
generalized neural network. It sweeps over ``N``, ``g`` and ``J0`` values, checks
for existing averages on disk, and computes only the missing ones. When the GPU
backend (``cupy``) is unavailable it falls back to ``numpy`` and parallelises the
number of averages across CPU cores.
"""

from __future__ import annotations

import argparse
import itertools
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from tqdm import tqdm

from kernels.SCS import SCSParameters, probe_output_graph
from kernels.SCS_TransCluster import (
    _compute_and_write,
    _compute_order_parameters,
    _make_build_output_path,
    _read_prev_partial,
)
from kernels.generic import resolve_backend, resolve_float_type, find_existing_averages
from lrgsglib.utils.basic.paths import remove_if_exists
from lrgsglib.utils.tools import ConditionalPartitioning


def _single_average_worker(task: tuple[dict, str, type, str, int | None, bool]) -> np.ndarray:
    """Compute one average sample and return [smax, smax2, gap, ediff]."""
    param_kwargs, backend, float_type, partition_rule, seed, use_edge_sign = task
    params = SCSParameters(seed=seed, **param_kwargs)
    partitioner = ConditionalPartitioning(partition_rule)
    gap, cluster_fraction, energy_diff = _compute_order_parameters(
        params, backend=backend, float_type=float_type, partitioner=partitioner,
        use_edge_sign_clustering=use_edge_sign
    )
    smax = cluster_fraction
    return np.array([smax, smax * smax, gap, energy_diff], dtype=float)


def _resolve_backend_with_fallback(backend: str) -> str:
    try:
        return resolve_backend(backend)
    except RuntimeError as exc:
        if backend.lower() == "cupy":
            print("Backend 'cupy' unavailable; falling back to 'numpy'.")
            return "numpy"
        raise exc


def _build_j0_values(j0_list: Sequence[float] | None, start: float, end: float, count: int) -> list[float]:
    if j0_list is not None and len(j0_list) > 0:
        return list(j0_list)
    return list(np.linspace(start, end, count))


def _accumulate_results(
    start_avg: int,
    accumulator_avg: np.ndarray,
    results: Iterable[np.ndarray],
) -> tuple[np.ndarray, int]:
    """Combine previous averages with new results."""
    totals = accumulator_avg * float(start_avg)
    n = start_avg
    for res in results:
        totals += res
        n += 1
    if n == 0:
        return accumulator_avg, 0
    return totals / float(n), n


def _process_point(
    param_kwargs: dict,
    *,
    backend: str,
    float_type: type,
    partition_rule: str,
    navg: int,
    out_suffix: str | None,
    nproc: int,
    combo_seed: int | None,
    verbose: bool,
    use_edge_sign_clustering: bool,
    pbar=None,
) -> None:
    params = SCSParameters(seed=None, **param_kwargs)
    probe = probe_output_graph(params)
    base_path = Path(probe.path_phtra)
    suffix_tokens = [
        "ordParam",
        probe.std_fname,
        f"J={params.J:.3g}",
        f"g={params.g:.3g}",
        f"J0={params.J0:.3g}",
    ]
    build_output_path = _make_build_output_path(base_path, suffix_tokens, out_suffix)

    pattern = "_".join(suffix_tokens) + "_na=*"
    if out_suffix:
        pattern += f"_{out_suffix}"
    pattern += ".txt"

    existing_navg = find_existing_averages(
        base_path,
        pattern,
        count_token="na",
        find_max=True,
        logger=None,
    )

    if existing_navg >= navg:
        if verbose and pbar:
            msg = f"[skip] N={params.N}, g={params.g}, J0={params.J0:.4g} already has {existing_navg} averages."
            pbar.set_postfix_str(msg, refresh=True)
        return

    accumulator_avg = np.zeros(4, dtype=float)
    if existing_navg > 0:
        prev_path = build_output_path(existing_navg)
        if prev_path.exists():
            try:
                accumulator_avg, prev_n = _read_prev_partial(prev_path)
                if prev_n != existing_navg and verbose and pbar:
                    msg = f"[warn] Expected navg={existing_navg} from filename but found {prev_n} inside {prev_path.name}"
                    pbar.set_postfix_str(msg, refresh=True)
            except Exception as exc:  # pragma: no cover - defensive
                if verbose and pbar:
                    msg = f"[warn] Could not read previous file {prev_path.name}: {exc}"
                    pbar.set_postfix_str(msg, refresh=True)
                accumulator_avg = np.zeros(4, dtype=float)
                existing_navg = 0
        else:
            if verbose and pbar:
                msg = f"[warn] Found navg={existing_navg} in filenames but file {prev_path.name} is missing; recomputing from scratch."
                pbar.set_postfix_str(msg, refresh=True)
            existing_navg = 0

    remaining = navg - existing_navg
    if remaining <= 0:
        return

    rng = np.random.default_rng(combo_seed)
    tasks = []
    for _ in range(remaining):
        sample_seed = int(rng.integers(0, 2**32))
        tasks.append((param_kwargs, backend, float_type, partition_rule, sample_seed, use_edge_sign_clustering))

    if nproc <= 1 or remaining == 1:
        iterator = map(_single_average_worker, tasks)
        if verbose and remaining > 1:
            iterator = tqdm(iterator, total=remaining, desc="averages (serial)", position=1, leave=False)
        results = list(iterator)
    else:
        with Pool(processes=nproc) as pool:
            iterator = pool.imap_unordered(_single_average_worker, tasks)
            if verbose:
                iterator = tqdm(iterator, total=remaining, desc="averages (pool)", position=1, leave=False)
            results = list(iterator)

    accumulator_avg, total_n = _accumulate_results(existing_navg, accumulator_avg, results)
    final_output_path = build_output_path(navg)
    _compute_and_write(final_output_path, accumulator_avg)

    if existing_navg > 0:
        remove_if_exists(build_output_path(existing_navg))

    if verbose and pbar:
        msg = f"[done] N={params.N}, g={params.g}, J0={params.J0:.4g} -> navg={total_n} written to {final_output_path.name}"
        pbar.set_postfix_str(msg, refresh=True)


def run_batch(
    *,
    N_list: Sequence[int],
    g_list: Sequence[float],
    j0_values: Sequence[float],
    gamma: float,
    J: float,
    diagonal: str,
    navg: int,
    backend: str,
    float_type: str,
    partition_rule: str,
    workdir: str,
    out_suffix: str | None,
    nproc: int | None,
    seed: int | None,
    verbose: bool,
    progress: bool = True,
    use_edge_sign_clustering: bool = False,
) -> None:
    backend_resolved = _resolve_backend_with_fallback(backend)
    float_type_resolved = resolve_float_type(float_type)
    if backend_resolved == "cupy":
        effective_nproc = 1
    elif nproc is None:
        effective_nproc = cpu_count()
    else:
        effective_nproc = max(1, nproc)

    combo_rng = np.random.default_rng(seed)

    combos = list(itertools.product(N_list, g_list, j0_values))

    pbar = None
    if progress:
        pbar = tqdm(total=len(combos), desc="sweep", position=0, leave=True)

    for N_val, g_val, j0_val in combos:
        param_kwargs = dict(
            N=int(N_val),
            J0=float(j0_val),
            gamma=float(gamma),
            J=float(J),
            g=float(g_val),
            diagonal=diagonal,
            workdir=workdir,
        )
        combo_seed = int(combo_rng.integers(0, 2**32))
        _process_point(
            param_kwargs,
            backend=backend_resolved,
            float_type=float_type_resolved,
            partition_rule=partition_rule,
            navg=navg,
            out_suffix=out_suffix,
            nproc=effective_nproc,
            combo_seed=combo_seed,
            verbose=verbose,
            use_edge_sign_clustering=use_edge_sign_clustering,
            pbar=pbar,
        )
        if pbar:
            pbar.update(1)

    if pbar:
        pbar.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep J0, g, and N for SCS order parameters (TSB-RB01 style).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-N",
        "--N",
        dest="N_list",
        type=int,
        nargs="+",
        default=[64],
        help="List of network sizes N.",
    )
    parser.add_argument(
        "--g",
        dest="g_list",
        type=float,
        nargs="+",
        default=[1.0],
        help="List of g coupling values.",
    )
    parser.add_argument(
        "--J0",
        dest="J0_list",
        type=float,
        nargs="+",
        default=None,
        help="Explicit J0 values (overrides linspace settings).",
    )
    parser.add_argument(
        "--J0-start",
        type=float,
        default=0.0,
        help="Start of the J0 linspace (used if --J0 is not provided).",
    )
    parser.add_argument(
        "--J0-end",
        type=float,
        default=3.0,
        help="End of the J0 linspace (used if --J0 is not provided).",
    )
    parser.add_argument(
        "--J0-count",
        type=int,
        default=61,
        help="Number of J0 samples in the linspace (used if --J0 is not provided).",
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=1.0,
        help="Gamma parameter for the SCS network.",
    )
    parser.add_argument(
        "--J",
        type=float,
        default=1.0,
        help="J coupling value.",
    )
    parser.add_argument(
        "--diagonal",
        type=str,
        default="-1.0",
        help="Diagonal specification (string or scalar, forwarded to SCSGeneralizedNN).",
    )
    parser.add_argument(
        "-na",
        "--navg",
        type=int,
        default=1000,
        help="Total number of averages to accumulate.",
    )
    parser.add_argument(
        "-wd",
        "--workdir",
        type=str,
        default="scs_ordParam_TSB-RB01",
        help="Working directory for storing results.",
    )
    parser.add_argument(
        "--partition-rule",
        type=str,
        default=">0",
        help="Partition rule for ConditionalPartitioning (e.g., '>0').",
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="cupy",
        help="Computation backend ('cupy' or 'numpy').",
    )
    parser.add_argument(
        "--float-type",
        type=str,
        default="float64",
        help="Floating point precision ('float32' or 'float64').",
    )
    parser.add_argument(
        "--out-suffix",
        type=str,
        default=None,
        help="Optional suffix appended to output files.",
    )
    parser.add_argument(
        "-np",
        "--nproc",
        type=int,
        default=None,
        help="Number of CPU processes to use when backend is numpy (default: all cores).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Base seed for reproducibility.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print progress information.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable the outer sweep progress bar.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    j0_values = _build_j0_values(args.J0_list, args.J0_start, args.J0_end, args.J0_count)

    run_batch(
        N_list=args.N_list,
        g_list=args.g_list,
        j0_values=j0_values,
        gamma=args.gamma,
        J=args.J,
        diagonal=args.diagonal,
        navg=args.navg,
        backend=args.backend,
        float_type=args.float_type,
        partition_rule=args.partition_rule,
        workdir=args.workdir,
        out_suffix=args.out_suffix,
        nproc=args.nproc,
        seed=args.seed,
        verbose=args.verbose,
        progress=not args.no_progress,
    )


if __name__ == "__main__":
    main()
