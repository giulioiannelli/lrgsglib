"""Kernel utilities for the 3D transient cluster simulations."""

from __future__ import annotations

import os
import pickle as pk
from collections import Counter
from typing import Callable, Iterable

import numpy as np

from lrgsglib import Lattice3D

from .L3D import make_transcluster_lattice
from .TransClusters import (
    build_geometry_selector,
    compose_output_path,
    load_pcluster_state,
    resolve_float_type,
    resolve_mode_extension,
)

__all__ = ["run_transcluster"]

# Debug flag - set to True to enable debug prints, False to disable
DEBUG = False


def _debug_print(*args, **kwargs):
    """Print debug messages only if DEBUG is enabled."""
    if DEBUG:
        print("[DEBUG]", *args, **kwargs)


def _probability_token(edge_weight: str, p_value: float, mu: float, sigma: float) -> str:
    if edge_weight == "flip":
        return f"p={p_value:.3g}"
    return f"mu={mu:.3g}_sigma={sigma:.3g}"


def file_path_maker(
    base_path: str,
    *,
    mode: str,
    p_value: float,
    navg: int | str,
    suffix: str,
    cell: str,
    ext: str,
    pdil: float,
    edge_weight: str,
    mu: float,
    sigma: float,
) -> str:
    """Compose the output file path for 3D simulations."""

    pdil_token = f"pdil={pdil:.3g}" if pdil else ""
    navg_token = f"na={navg}" if navg not in {"", None} else ""
    probability = _probability_token(edge_weight, p_value, mu, sigma)
    return compose_output_path(
        base_path,
        mode,
        pdil_token,
        probability,
        cell,
        navg_token,
        suffix,
        ext=ext,
    )


def _instantiate_lattice(args, *, init_nw_dict: bool, with_positions: bool = False) -> Lattice3D:
    return make_transcluster_lattice(
        args,
        init_nw_dict=init_nw_dict,
        with_positions=with_positions,
        only_const_mode=False,
    )


def run_pcluster_mode(
    args,
    geometry_func: Callable[[Lattice3D], Iterable],
    base_path: str,
    ext: str,
    save_frequency: int,
    typf,
) -> None:
    if args.edge_weight != "flip":
        raise ValueError("pCluster mode requires 'flip' edge weights in 3D simulations.")

    base_pattern = file_path_maker(
        base_path,
        mode=args.mode,
        p_value=args.p,
        navg="",
        suffix="",
        cell=args.cell_type,
        ext="",
        pdil=args.pdil,
        edge_weight=args.edge_weight,
        mu=args.mu,
        sigma=args.sigma,
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
            pdil=args.pdil,
            edge_weight=args.edge_weight,
            mu=args.mu,
            sigma=args.sigma,
        )

    averages_completed = n_avg_done
    while averages_completed < args.number_of_averages:
        lattice = _instantiate_lattice(args, init_nw_dict=True)
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
            pdil=args.pdil,
            edge_weight=args.edge_weight,
            mu=args.mu,
            sigma=args.sigma,
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
    geometry_func: Callable[[Lattice3D], Iterable],
    base_path: str,
    ext: str,
    save_frequency: int,
    typf,
) -> None:
    pinf: list[float] = []
    pinf_sq: list[float] = []
    neglinks = 0
    last_written: str | None = None
    last_data: list[float] | None = None

    def persist(avg_count: int, data: list[float]) -> None:
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
            pdil=args.pdil,
            edge_weight=args.edge_weight,
            mu=args.mu,
            sigma=args.sigma,
        )
        with open(filename, "wb") as handle:
            np.savetxt(handle, np.atleast_2d(data), fmt="%.7g")
        last_written = filename

    for avg_idx in range(args.number_of_averages):
        lattice = _instantiate_lattice(args, init_nw_dict=True)
        if args.edge_weight == "flip":
            lattice.flip_sel_edges(geometry_func(lattice))
        else:
            lattice.set_edges_random_normal(args.mu, args.sigma)

        neglinks += lattice.Ne_n
        lattice.compute_k_eigvV(typf=typf)

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
            pinf_var = max(pinf_sq_mean - pinf_mean**2, 0.0)
            pinf_std = float(np.sqrt(pinf_var))
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

    preview = make_transcluster_lattice(
        args,
        init_nw_dict=False,
        with_positions=False,
        only_const_mode=True,
    )
    base_path = {
        "pCluster": getattr(preview, "path_lrgsg"),
        "ordParam": getattr(preview, "path_phtra"),
    }[args.mode]

    filename = file_path_maker(
        base_path,
        mode=args.mode,
        p_value=args.p,
        navg=navg,
        suffix=args.out_suffix,
        cell=args.cell_type,
        ext=ext,
        pdil=args.pdil,
        edge_weight=args.edge_weight,
        mu=args.mu,
        sigma=args.sigma,
    )
    if os.path.exists(filename):
        raise SystemExit(f"File {os.path.split(filename)[1]} already exists.")

    if args.mode == "pCluster":
        run_pcluster_mode(args, geometry_func, base_path, ext, save_frequency, typf)
    else:
        run_ordparam_mode(args, geometry_func, base_path, ext, save_frequency, typf)
