"""Computational kernels for the 2D transient cluster analysis program."""

from __future__ import annotations

import glob
import os
import pickle as pk
import re
from collections import Counter
from pathlib import Path
from typing import Any, Callable

import numpy as np

from lrgsglib.config.const import PKL, TXT
from lrgsglib.core import Lattice2D
from lrgsglib.utils.basic.strings import get_first_int_in_str
from .L2D import initialize_l2d_dict_args


def resolve_float_type(float_type: str) -> type[np.floating]:
    """Convert a CLI string argument into a NumPy floating type."""

    mapping = {
        "float32": np.float32,
        "float64": np.float64,
    }
    try:
        return mapping[float_type]
    except KeyError as exc:
        raise ValueError("Invalid float type specified") from exc


def resolve_mode_extension(mode: str) -> str:
    """Return the output file extension associated with a simulation mode."""

    match mode:
        case "pCluster":
            return PKL
        case "ordParam":
            return TXT
        case _:
            raise ValueError("Invalid mode specified")


def build_geometry_func(cell: str) -> Callable[[Lattice2D], list[int]]:
    """Return the geometry selector associated with the required cell type."""

    match cell:
        case "rand" | "randZERR" | "randXERR":

            def geometry_func(lattice: Lattice2D) -> list[int]:
                return lattice.nwDict[cell]["G"]

        case _ if cell.startswith("ball"):
            radius = get_first_int_in_str(cell)

            def geometry_func(lattice: Lattice2D) -> list[int]:
                return lattice.nwDict.get_links_rball(radius)

        case _:
            raise ValueError("Invalid cell specified")

    return geometry_func


def build_file_path(
    base_path: str | os.PathLike[str],
    *,
    mode: str,
    p_value: float,
    averages: int | str,
    cell_type: str,
    prewiring: float,
    suffix: str,
    extension: str | None = None,
) -> str:
    """Construct the full path to an output artefact for the requested run."""

    ext = extension if extension is not None else resolve_mode_extension(mode)
    prewiring_str = f"prew={prewiring:.3g}" if prewiring != 0.0 else ""
    name_tokens = filter(
        None,
        [mode, f"p={p_value:.3g}", cell_type, f"na={averages}", prewiring_str, suffix],
    )
    filename = "_".join(name_tokens) + ext
    return os.path.join(base_path, filename)


def ensure_output_path_available(filepath: str) -> None:
    """Abort execution if the destination file already exists."""

    if os.path.exists(filepath):
        raise SystemExit(f"File {Path(filepath).name} already exists.")


def build_lattice_kwargs(
    args: Any,
    *,
    with_positions: bool = False,
    require_nw_dict: bool = False,
) -> dict[str, Any]:
    """Return keyword arguments for ``Lattice2D`` consistent with parser defaults."""

    kwargs = initialize_l2d_dict_args(args)
    kwargs["sgpathn"] = args.workdir
    if require_nw_dict:
        kwargs["init_nw_dict"] = True
    if with_positions:
        kwargs["with_positions"] = True
    if args.prew:
        kwargs["prew"] = args.prew
    return kwargs


def create_probe_lattice(args) -> Lattice2D:
    """Instantiate a lattice to determine output directories."""

    return Lattice2D(**build_lattice_kwargs(args))


def determine_output_paths(lattice: Lattice2D) -> dict[str, str]:
    """Return the mapping between simulation modes and their output folders."""

    return {
        "pCluster": lattice.path_lrgsg,
        "ordParam": lattice.path_phtra.name,
    }


def pcluster_mode(
    args,
    *,
    geometry_func: Callable[[Lattice2D], list[int]],
    base_path: str,
    dtype: type[np.floating],
    save_frequency: int,
) -> None:
    """Execute the pCluster mode simulation."""

    merged_dict: Counter[int] = Counter()
    n_avg_done = 0
    file_pattern = build_file_path(
        base_path,
        mode=args.mode,
        p_value=args.p,
        averages="",
        cell_type=args.cell_type,
        prewiring=args.prew,
        suffix="",
        extension="",
    )

    try:
        existing_file = glob.glob(f"{file_pattern}*")[0]
        with open(existing_file, "rb") as file:
            merged_dict = pk.load(file)

        avg_index = -2 if args.out_suffix else -1
        n_avg_done = os.path.splitext(existing_file.split("_")[avg_index])[0]
        n_avg_done = int(re.search(r"\d+", n_avg_done).group())
        previous_filename = existing_file
    except Exception:  # noqa: BLE001 - Preserve legacy broad exception behaviour
        previous_filename = build_file_path(
            base_path,
            mode=args.mode,
            p_value=args.p,
            averages=0,
            cell_type=args.cell_type,
            prewiring=args.prew,
            suffix=args.out_suffix,
        )

    averages_needed = args.number_of_averages - n_avg_done
    batches = (averages_needed // save_frequency) + bool(averages_needed % save_frequency)

    for current_period in range(batches):
        batch_size = min(averages_needed - current_period * save_frequency, save_frequency)
        for _ in range(batch_size):
            lattice = Lattice2D(
                **build_lattice_kwargs(args, require_nw_dict=True),
            )
            lattice.flip_sel_edges(geometry_func(lattice))
            lattice.compute_k_eigvV(typf=dtype)
            dist_dict = lattice.get_cluster_distribution()
            merged_dict += dist_dict

        averages_completed = n_avg_done + (current_period + 1) * save_frequency
        new_filename = build_file_path(
            base_path,
            mode=args.mode,
            p_value=args.p,
            averages=averages_completed,
            cell_type=args.cell_type,
            prewiring=args.prew,
            suffix=args.out_suffix,
        )
        try:
            os.rename(previous_filename, new_filename)
        except (FileNotFoundError, OSError):
            pass

        with open(new_filename, "wb") as file:
            pk.dump(merged_dict, file)
        previous_filename = new_filename


def order_parameter_mode(
    args,
    *,
    geometry_func: Callable[[Lattice2D], list[int]],
    base_path: str,
    dtype: type[np.floating],
    save_frequency: int,
) -> None:
    """Execute the order-parameter accumulation workflow."""

    pinf_values: list[float] = []
    pinf_variance: list[float] = []
    negative_links = 0
    averages_completed = 0

    for _ in range(args.number_of_averages):
        lattice = Lattice2D(
            **build_lattice_kwargs(
                args,
                with_positions=True,
                require_nw_dict=True,
            ),
        )
        lattice.flip_sel_edges(geometry_func(lattice))
        lattice.compute_k_eigvV(typf=dtype)

        averages_completed += 1
        negative_links += lattice.Ne_n

        lattice.load_eigV_on_graph(binarize=True)
        lattice.make_clustersYN("eigV0", +1)

        if len(lattice.clustersY) > 1:
            smax2 = len(max(lattice.clustersY, key=len)) / (1.0 * lattice.N)
        else:
            smax2 = 1.0

        pinf_values.append(smax2)
        pinf_variance.append(smax2)

        data = [
            averages_completed,
            lattice.pflip,
            negative_links / averages_completed,
            np.mean(pinf_values),
            np.mean(pinf_variance),
            np.std(pinf_values),
        ]

        if averages_completed % save_frequency == 0:
            try:
                previous_file = build_file_path(
                    base_path,
                    mode=args.mode,
                    p_value=args.p,
                    averages=averages_completed - save_frequency,
                    cell_type=args.cell_type,
                    prewiring=args.prew,
                    suffix=args.out_suffix,
                )
                os.remove(previous_file)
            except OSError:
                pass

            filename = build_file_path(
                base_path,
                mode=args.mode,
                p_value=args.p,
                averages=averages_completed,
                cell_type=args.cell_type,
                prewiring=args.prew,
                suffix=args.out_suffix,
            )
            with open(filename, "wb") as file:
                np.savetxt(file, np.atleast_2d(data), fmt="%.7g")

    remains = args.number_of_averages % save_frequency
    os.rename(
        build_file_path(
            base_path,
            mode=args.mode,
            p_value=args.p,
            averages=args.number_of_averages - remains,
            cell_type=args.cell_type,
            prewiring=args.prew,
            suffix=args.out_suffix,
        ),
        build_file_path(
            base_path,
            mode=args.mode,
            p_value=args.p,
            averages=args.number_of_averages,
            cell_type=args.cell_type,
            prewiring=args.prew,
            suffix=args.out_suffix,
        ),
    )


def execute_transcluster(args) -> None:
    """Run the appropriate transient-cluster routine based on CLI arguments."""

    computed_frequency = args.number_of_averages // 20
    save_frequency = args.save_frequency or max(1, computed_frequency)
    dtype = resolve_float_type(args.float_type)
    geometry_func = build_geometry_func(args.cell_type)

    probe_lattice = create_probe_lattice(args)
    mode_paths = determine_output_paths(probe_lattice)

    base_path = mode_paths.get(args.mode)
    if base_path is None:
        raise ValueError("Invalid mode specified")

    output_filepath = build_file_path(
        base_path,
        mode=args.mode,
        p_value=args.p,
        averages=args.number_of_averages,
        cell_type=args.cell_type,
        prewiring=args.prew,
        suffix=args.out_suffix,
    )
    ensure_output_path_available(output_filepath)

    if args.mode == "pCluster":
        pcluster_mode(
            args,
            geometry_func=geometry_func,
            base_path=base_path,
            dtype=dtype,
            save_frequency=save_frequency,
        )
    elif args.mode == "ordParam":
        order_parameter_mode(
            args,
            geometry_func=geometry_func,
            base_path=base_path,
            dtype=dtype,
            save_frequency=save_frequency,
        )
    else:
        raise ValueError("Invalid mode specified")
