"""Kernel utilities for the generalized SCS transient cluster simulations."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from lrgsglib import SCSGeneralizedNN
from lrgsglib.utils.tools import ConditionalPartitioning

from .SCS import (
    SCSParameters,
    build_scs_kwargs,
    probe_output_graph,
    resolve_backend,
    resolve_float_type,
)

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
) -> tuple[float, float, float, float]:
    kwargs = build_scs_kwargs(params, J0_value=j0, make_dir_tree=False, seed=seed)
    scs = SCSGeneralizedNN(**kwargs)
    scs.compute_laplacian_spectrum_weigV(backend=backend, typf=float_type)
    scs.load_eigV_on_graph(which=0, binarize=False)
    scs.make_eigVclustersYN(val=partitioner, which=0, binarize=True)

    gap = (scs.eigv[1] - scs.eigv[0]) / scs.eigv[-1]
    cluster_fraction = _largest_cluster_fraction(scs)
    energy0 = scs.get_sksph_energy_eigV(0)
    energy1 = scs.get_sksph_energy_eigV(1)
    return gap, cluster_fraction, cluster_fraction * cluster_fraction, energy0 - energy1


def _write_results(path: Path, data: Iterable[Iterable[float]]) -> None:
    array = np.asarray(list(data), dtype=float)
    if array.size == 0:
        return
    np.savetxt(path, array, fmt="%.7g")


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
        f"navg={args.number_of_averages}",
    ]
    if args.out_suffix:
        suffix_tokens.append(args.out_suffix)
    filename = "_".join(suffix_tokens) + ".txt"
    output_path = Path(base_dir) / filename

    del probe

    if output_path.exists():
        if getattr(args, "remove_files", False):
            output_path.unlink()
        else:
            raise SystemExit(f"File {output_path.name} already exists.")

    results: list[list[float]] = []
    rng = np.random.default_rng(args.seed) if args.seed is not None else None

    accumulator = np.zeros(4, dtype=float)
    for avg_idx in range(args.number_of_averages):
        sample_seed = int(rng.integers(0, 2**32)) if rng is not None else None
        gap, smax, smax2, energy_diff = _compute_order_parameters(
            params,
            j0=j0_value,
            backend=backend,
            float_type=float_type,
            partitioner=partitioner,
            seed=sample_seed,
        )
        accumulator += np.array([gap, smax, smax2, energy_diff])
        if args.save_frequency and args.save_frequency > 0:
            if (avg_idx + 1) % args.save_frequency == 0:
                interim_average = accumulator / float(avg_idx + 1)
                _write_results(output_path, [[j0_value, *interim_average.tolist()]])

    accumulator /= args.number_of_averages
    results.append([j0_value, *accumulator.tolist()])
    _write_results(output_path, results)
