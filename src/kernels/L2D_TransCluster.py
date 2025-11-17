"""Kernel utilities for the 2D transient cluster simulations."""

from __future__ import annotations

import os
import pickle as pk
from collections import Counter
from typing import Callable, Iterable

import numpy as np

from lrgsglib import Lattice2D

from .TransClusters import (
    build_geometry_selector,
    compose_output_path,
    load_pcluster_state,
    resolve_mode_extension,
)
from .generic import resolve_float_type

__all__ = ["run_transcluster"]


def file_path_maker(
    base_path: str,
    *,
    mode: str,
    p_value: float,
    navg: int | str,
    suffix: str,
    cell: str,
    ext: str,
    prew: float,
) -> str:
    """Compose the output file path following the legacy naming scheme."""

    prew_str = f"prew={prew:.3g}" if prew != 0.0 else ""
    navg_token = f"na={navg}" if navg not in {"", None} else ""
    return compose_output_path(
        base_path,
        mode,
        f"p={p_value:.3g}",
        cell,
        navg_token,
        prew_str,
        suffix,
        ext=ext,
    )


def run_pcluster_mode(
    args,
    geometry_func: Callable[[Lattice2D], Iterable],
    base_path: str,
    ext: str,
    save_frequency: int,
    typf,
) -> None:
    """Execute the ``pCluster`` averaging routine."""

    base_pattern = file_path_maker(
        base_path,
        mode=args.mode,
        p_value=args.p,
        navg="",
        suffix="",
        cell=args.cell_type,
        ext="",
        prew=args.prew,
    )
    merged_dict, n_avg_done, fname_old = load_pcluster_state(
        base_pattern, bool(args.out_suffix)
    )
    merged_dict = Counter(merged_dict)

    if fname_old is None:
        fname_old = file_path_maker(
            base_path,
            mode=args.mode,
            p_value=args.p,
            navg=0,
            suffix=args.out_suffix,
            cell=args.cell_type,
            ext=ext,
            prew=args.prew,
        )

    averages_completed = n_avg_done
    while averages_completed < args.number_of_averages:
        lattice = Lattice2D(
            args.L,
            pflip=args.p,
            geo=args.geometry,
            init_nw_dict=True,
        )
        lattice.flip_sel_edges(geometry_func(lattice))
        lattice.compute_k_eigvV(typf=typf)
        merged_dict += lattice.get_cluster_distribution()

        averages_completed += 1
        should_store = (
            (averages_completed - n_avg_done) % save_frequency == 0
            or averages_completed == args.number_of_averages
        )
        if not should_store:
            continue

        fname_new = file_path_maker(
            base_path,
            mode=args.mode,
            p_value=args.p,
            navg=averages_completed,
            suffix=args.out_suffix,
            cell=args.cell_type,
            ext=ext,
            prew=args.prew,
        )
        if fname_old and fname_old != fname_new:
            try:
                os.rename(fname_old, fname_new)
            except (FileNotFoundError, OSError):
                pass
        with open(fname_new, "wb") as handle:
            pk.dump(merged_dict, handle)
        fname_old = fname_new


def run_ordparam_mode(
    args,
    geometry_func: Callable[[Lattice2D], Iterable],
    base_path: str,
    ext: str,
    save_frequency: int,
    typf,
) -> None:
    """Execute the ``ordParam`` averaging routine."""

    pinf: list[float] = []
    pinf_sq: list[float] = []
    neglinks = 0
    last_written: str | None = None
    last_data: list[float] | None = None

    def persist(avg_count: int, data: list[float]) -> str:
        nonlocal last_written
        if last_written:
            try:
                os.remove(last_written)
            except OSError:
                pass

        filename = file_path_maker(
            base_path,
            mode=args.mode,
            p_value=args.p,
            navg=avg_count,
            suffix=args.out_suffix,
            cell=args.cell_type,
            ext=ext,
            prew=args.prew,
        )
        with open(filename, "wb") as handle:
            np.savetxt(handle, np.atleast_2d(data), fmt="%.7g")
        last_written = filename
        return filename

    for avg_idx in range(args.number_of_averages):
        lattice = Lattice2D(
            args.L,
            pflip=args.p,
            geo=args.geometry,
            with_positions=True,
            init_nw_dict=True,
            sgpathn=args.workdir,
        )
        lattice.flip_sel_edges(geometry_func(lattice))
        lattice.compute_k_eigvV(typf=typf)

        neglinks += lattice.Ne_n

        try:
            lattice.compute_pinf()
        except IndexError:
            continue

        pinf_val = float(lattice.Pinf)
        pinf.append(pinf_val)
        pinf_sq.append(pinf_val * pinf_val)

        avg_count = avg_idx + 1
        if pinf:
            pinf_mean = float(np.mean(pinf))
            pinf_sq_mean = float(np.mean(pinf_sq))
            pinf_std = float(np.sqrt(max(pinf_sq_mean - pinf_mean**2, 0.0)))
        else:
            pinf_mean = pinf_sq_mean = pinf_std = 0.0

        data = [
            avg_count,
            lattice.pflip,
            neglinks / avg_count,
            pinf_mean,
            pinf_sq_mean,
            pinf_std,
        ]

        last_data = data

        if save_frequency and avg_count % save_frequency == 0:
            persist(avg_count, data)

    if args.number_of_averages > 0 and last_data is not None:
        persist(args.number_of_averages, last_data)


def run_transcluster(args) -> None:
    """Driver used by the CLI entry-point."""

    if args.number_of_averages <= 0:
        raise ValueError("Number of averages must be positive.")

    navg = args.number_of_averages
    save_frequency = args.save_frequency if args.save_frequency else max(1, navg // 20)

    typf = resolve_float_type(args.float_type)
    ext = resolve_mode_extension(args.mode)
    geometry_func = build_geometry_selector(args.cell_type)

    test_lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry, sgpathn=args.workdir)
    base_path = {
        "pCluster": getattr(test_lattice, "path_lrgsg"),
        "ordParam": getattr(test_lattice, "path_phtra"),
    }[args.mode]

    filename = file_path_maker(
        base_path,
        mode=args.mode,
        p_value=args.p,
        navg=navg,
        suffix=args.out_suffix,
        cell=args.cell_type,
        ext=ext,
        prew=args.prew,
    )
    if os.path.exists(filename):
        raise SystemExit(f"File {os.path.split(filename)[1]} already exists.")

    if args.mode == "pCluster":
        run_pcluster_mode(args, geometry_func, base_path, ext, save_frequency, typf)
    else:
        run_ordparam_mode(args, geometry_func, base_path, ext, save_frequency, typf)

