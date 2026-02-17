"""Kernel for Ising dynamics on 3D lattices.

Combines topology helpers from :mod:`kernels.L3D` with dynamics helpers
from :mod:`kernels.IsingDynamics`.  Engine-agnostic: the graph engine
(NX vs GT) is handled entirely by :func:`~kernels.L3D.prepare_lattice`.
"""

from __future__ import annotations

from typing import Any

from lrgsglib import *
from .L3D import get_l3d_out_suffix, prepare_lattice
from .IsingDynamics import (
    initialize_ising_dict_args,
    run_ising_dynamics,
    clean_up_files,
)

__all__ = ["run_simulation"]


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
            if args.remove_files:
                clean_up_files(isdy, lattice)
