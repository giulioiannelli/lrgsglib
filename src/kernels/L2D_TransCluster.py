"""Kernel utilities for the 2D transient cluster simulations."""

from __future__ import annotations

import glob
import os
import pickle as pk
import re
from collections import Counter
from typing import Callable, Tuple

import numpy as np

from lrgsglib import Lattice2D
from lrgsglib.config.const import PKL, TXT
from lrgsglib.utils.basic.strings import get_first_int_in_str

__all__ = ["run_transcluster"]


def resolve_float_type(float_type: str) -> type:
    """Map a float type string to the corresponding NumPy dtype."""
    match float_type:
        case "float32":
            return np.float32
        case "float64":
            return np.float64
        case _:
            raise ValueError("Invalid float type specified")


def resolve_mode_extension(mode: str) -> str:
    """Return the output file extension associated with a computation mode."""
    match mode:
        case "pCluster":
            return PKL
        case "ordParam":
            return TXT
        case _:
            raise ValueError("Invalid mode specified")


def get_geometry_func(cell: str) -> Callable[[Lattice2D], np.ndarray]:
    """Return the geometry selection function based on the defect type."""
    match cell:
        case "rand" | "randZERR" | "randXERR":
            def geometry_func(lattice: Lattice2D):
                return lattice.nwDict[cell]["G"]
        case _ if cell.startswith("ball"):
            radius = get_first_int_in_str(cell)

            def geometry_func(lattice: Lattice2D):
                return lattice.nwDict.get_links_rball(radius)
        case _:
            raise ValueError("Invalid cell specified")

    return geometry_func


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
    str_name = "_".join(
        filter(
            None,
            [mode, f"p={p_value:.3g}", cell, f"na={navg}", prew_str, suffix],
        )
    )
    return os.path.join(base_path, str_name) + ext


def load_existing_pcluster_state(
    base_pattern: str,
    has_suffix: bool,
) -> Tuple[Counter, int, str | None]:
    """Load the latest persisted cluster distribution if present."""
    try:
        fname_exists = glob.glob(f"{base_pattern}*")[0]
    except IndexError:
        return Counter(), 0, None

    with open(fname_exists, "rb") as handle:
        merged_dict = pk.load(handle)

    avg_idx = -2 if has_suffix else -1
    avg_token = os.path.splitext(fname_exists.split("_")[avg_idx])[0]
    n_avg_done = int(re.search(r"\d+", avg_token).group())
    return merged_dict, n_avg_done, fname_exists


def run_pcluster_mode(
    args,
    geometry_func: Callable[[Lattice2D], np.ndarray],
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
    merged_dict, n_avg_done, fname_old = load_existing_pcluster_state(
        base_pattern, bool(args.out_suffix)
    )
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
    n_avg_need = args.number_of_averages - n_avg_done

    for current_period in range((n_avg_need // save_frequency) + bool(n_avg_need % save_frequency)):
        batch_size = min(
            n_avg_need - current_period * save_frequency,
            save_frequency,
        )
        for _ in range(batch_size):
            lattice = Lattice2D(
                args.L,
                pflip=args.p,
                geo=args.geometry,
                init_nw_dict=True,
            )
            lattice.flip_sel_edges(geometry_func(lattice))
            lattice.compute_k_eigvV(typf=typf)
            merged_dict += lattice.get_cluster_distribution()

        navg_curr = n_avg_done + (current_period + 1) * save_frequency
        fname_new = file_path_maker(
            base_path,
            mode=args.mode,
            p_value=args.p,
            navg=navg_curr,
            suffix=args.out_suffix,
            cell=args.cell_type,
            ext=ext,
            prew=args.prew,
        )
        try:
            os.rename(fname_old, fname_new)
        except (FileNotFoundError, OSError):
            pass
        with open(fname_new, "wb") as handle:
            pk.dump(merged_dict, handle)
        fname_old = fname_new


def run_ordparam_mode(
    args,
    geometry_func: Callable[[Lattice2D], np.ndarray],
    base_path: str,
    ext: str,
    save_frequency: int,
    typf,
) -> None:
    """Execute the ``ordParam`` averaging routine."""
    pinf: list[float] = []
    pinf_var: list[float] = []
    neglinks = 0
    avg_count = 0

    for _ in range(args.number_of_averages):
        lattice = Lattice2D(
            args.L,
            pflip=args.p,
            geo=args.geometry,
            with_positions=True,
            init_nw_dict=True,
            sgpathn=args.outdir,
        )
        lattice.flip_sel_edges(geometry_func(lattice))
        lattice.compute_k_eigvV(typf=typf)

        avg_count += 1
        neglinks += lattice.Ne_n

        lattice.load_eigV_on_graph(binarize=True)
        lattice.make_clustersYN("eigV0", +1)
        if len(lattice.clustersY) > 1:
            smax2 = len(max(lattice.clustersY, key=len)) / (1.0 * lattice.N)
        else:
            smax2 = 1.0

        pinf.append(smax2)
        pinf_var.append(smax2)

        data = [
            avg_count,
            lattice.pflip,
            neglinks / avg_count,
            np.mean(pinf),
            np.mean(pinf_var),
            np.std(pinf),
        ]

        if avg_count % save_frequency == 0:
            try:
                old_file = file_path_maker(
                    base_path,
                    mode=args.mode,
                    p_value=args.p,
                    navg=avg_count - save_frequency,
                    suffix=args.out_suffix,
                    cell=args.cell_type,
                    ext=ext,
                    prew=args.prew,
                )
                os.remove(old_file)
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

    remains = args.number_of_averages % save_frequency
    os.rename(
        file_path_maker(
            base_path,
            mode=args.mode,
            p_value=args.p,
            navg=args.number_of_averages - remains,
            suffix=args.out_suffix,
            cell=args.cell_type,
            ext=ext,
            prew=args.prew,
        ),
        file_path_maker(
            base_path,
            mode=args.mode,
            p_value=args.p,
            navg=args.number_of_averages,
            suffix=args.out_suffix,
            cell=args.cell_type,
            ext=ext,
            prew=args.prew,
        ),
    )


def run_transcluster(args) -> None:
    """Driver used by the CLI entry-point."""
    navg = args.number_of_averages
    save_frequency = args.save_frequency if args.save_frequency else navg // 20

    typf = resolve_float_type(args.float_type)
    ext = resolve_mode_extension(args.mode)
    geometry_func = get_geometry_func(args.cell_type)

    test_lattice = Lattice2D(args.L, pflip=args.p, geo=args.geometry, sgpathn=args.outdir)
    base_path = {
        "pCluster": test_lattice.path_lrgsg,
        "ordParam": test_lattice.path_phtra.name,
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
        exit(f"File {os.path.split(filename)[1]} already exists.")

    if args.mode == "pCluster":
        run_pcluster_mode(args, geometry_func, base_path, ext, save_frequency, typf)
    else:
        run_ordparam_mode(args, geometry_func, base_path, ext, save_frequency, typf)
