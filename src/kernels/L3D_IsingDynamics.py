"""Kernel for Ising dynamics on 3D lattices.

Combines topology helpers from :mod:`kernels.L3D` with dynamics helpers
from :mod:`kernels.IsingDynamics`.  Engine-agnostic: the graph engine
(NX vs GT) is handled entirely by :func:`~kernels.L3D.prepare_lattice`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from lrgsglib import *
from .L3D import get_l3d_out_suffix, prepare_lattice
from .IsingDynamics import (
    initialize_ising_dict_args,
    run_ising_dynamics,
    clean_up_files,
)

__all__ = ["run_simulation"]


def _save_ising_results(isdy, lattice, args, quench_idx: int) -> Path:
    """Save energy/magnetization trajectories to an NPZ file.

    Parameters
    ----------
    isdy : IsingDynamics
        Executed IsingDynamics instance with results.
    lattice : SignedGraph
        The lattice used for the simulation.
    args : argparse.Namespace
        Parsed CLI arguments.
    quench_idx : int
        Index of the current disorder realization.

    Returns
    -------
    Path
        Path to the saved NPZ file.
    """
    save_dir = Path(lattice.path_sgdata) / "ising_results"
    save_dir.mkdir(parents=True, exist_ok=True)

    rl = isdy.runlang
    side = getattr(lattice, "side1", 0)
    pflip = getattr(lattice, "pflip", 0.0)
    fname = f"ising_{rl}_L{side}_p{pflip:.3g}_q{quench_idx:03d}.npz"

    data: dict[str, Any] = {
        "runlang": rl,
        "side": side,
        "pflip": pflip,
        "N": isdy.N,
    }
    # Standard trajectories
    if hasattr(isdy, "ene") and isdy.ene is not None:
        data["ene"] = np.asarray(isdy.ene)
    if hasattr(isdy, "magn") and isdy.magn is not None:
        data["magn"] = np.asarray(isdy.magn)
    # SA trajectories
    if hasattr(isdy, "sa_energy") and isdy.sa_energy is not None:
        data["sa_energy"] = np.asarray(isdy.sa_energy)
        data["sa_magn"] = np.asarray(isdy.sa_magn)
        data["sa_temps"] = np.asarray(isdy.sa_temps)
    # Topo_met results
    if hasattr(isdy, "topo_met_best_energy"):
        data["topo_met_best_energy"] = isdy.topo_met_best_energy
        if hasattr(isdy, "topo_met_best_spins"):
            data["topo_met_best_spins"] = isdy.topo_met_best_spins
        if hasattr(isdy, "topo_met_coeffs"):
            data["topo_met_coeffs"] = isdy.topo_met_coeffs
    # Topo_fca results
    if hasattr(isdy, "topo_fca_field") and isdy.topo_fca_field is not None:
        data["topo_fca_field"] = isdy.topo_fca_field

    out_path = save_dir / fname
    np.savez_compressed(out_path, **data)
    return out_path


def run_simulation(args: Any) -> None:
    """Run L3D Ising dynamics with quenched + thermal averaging.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments containing lattice, Ising, and averaging
        parameters.
    """
    n_quench = getattr(args, "number_of_averages", 1)
    n_thermal = getattr(args, "n_thermal", 1)
    out_suffix = get_l3d_out_suffix(args)

    for _q in range(n_quench):
        lattice = prepare_lattice(args)

        # Optional: eigenvector-based init for ground_state_k conditions
        # Requires compute_k_eigvV / load_eigV_on_graph on the graph class.
        ic = getattr(args, "init_cond", "rand")
        if ic.startswith("ground_state") or ic.startswith("gs"):
            k = int(ic.rsplit("_", 1)[-1])
            if hasattr(lattice, "compute_k_eigvV") and hasattr(
                lattice, "load_eigV_on_graph"
            ):
                lattice.compute_k_eigvV(k + 1)
                lattice.load_eigV_on_graph(k, binarize=True)
            else:
                import warnings

                warnings.warn(
                    f"ground_state init requested but {type(lattice).__name__} "
                    "does not support eigenvector loading. Using random init.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                args.init_cond = "rand"

        for _t in range(n_thermal):
            isdy = run_ising_dynamics(args, lattice, out_suffix=out_suffix)
            if getattr(args, 'save_results', False):
                _save_ising_results(
                    isdy, lattice, args, _q * n_thermal + _t
                )
            if args.remove_files:
                clean_up_files(isdy, lattice)
