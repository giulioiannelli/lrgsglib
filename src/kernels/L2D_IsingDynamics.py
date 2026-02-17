from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from lrgsglib import *
from .L2D import *
from .IsingDynamics import *
from parsers.shared import get_graph_engine


def _save_ising_results(isdy, lattice, args, quench_idx: int) -> Path:
    """Save energy/magnetization trajectories to an NPZ file (L2D).

    Parameters
    ----------
    isdy : IsingDynamics
        Executed IsingDynamics instance with results.
    lattice : SignedGraph
        The 2D lattice used for the simulation.
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
    geo = getattr(lattice, "geo", "sqr")
    fname = (
        f"ising_{rl}_L{side}_{geo}_p{pflip:.3g}_q{quench_idx:03d}.npz"
    )

    data: dict[str, Any] = {
        "runlang": rl,
        "side": side,
        "pflip": pflip,
        "geo": geo,
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


def run_simulation(args):
    engine = get_graph_engine(args)

    if engine == "gt":
        _run_simulation_gt(args)
    else:
        _run_simulation_nx(args)


def _run_simulation_gt(args):
    """GT path: skip eigvec/cluster pipeline, use direct Ising dynamics."""
    for _q in range(args.number_of_averages):
        lattice = prepare_lattice(args)
        isdy = run_ising_dynamics(args, lattice)
        if getattr(args, 'save_results', False):
            _save_ising_results(isdy, lattice, args, _q)
        if args.remove_files:
            clean_up_files(isdy, lattice)


def _run_simulation_nx(args):
    """NX path: eigenvector/cluster-based init for ground_state, direct for others."""
    ic_gs = args.init_cond.startswith('ground_state') or args.init_cond.startswith('gs')

    if not ic_gs:
        # Direct path: no eigenvector/cluster pipeline needed
        for _q in range(args.number_of_averages):
            lattice = Lattice2D(**initialize_l2d_dict_args(args))
            if lattice.init_nw_dict:
                lattice.flip_sel_edges(lattice.nwDict[args.cell_type]['G'])
            else:
                lattice.flip_random_fract_edges()
            isdy = run_ising_dynamics(args, lattice)
            if getattr(args, 'save_results', False):
                _save_ising_results(isdy, lattice, args, _q)
            if args.remove_files:
                clean_up_files(isdy, lattice)
        return

    # Ground-state path: eigenvector/cluster-based initialization
    number = int(args.init_cond.split('_')[-1])
    val = ConditionalPartitioning(args.val)
    for _q in range(args.number_of_averages):
        lattice = Lattice2D(**initialize_l2d_dict_args(args))
        if lattice.init_nw_dict:
            lattice.flip_sel_edges(lattice.nwDict[args.cell_type]['G'])
        else:
            lattice.flip_random_fract_edges()
        lattice.compute_k_eigvV(number+1)
        lattice.load_eigV_on_graph(number, binarize=True)
        lattice.make_clustersYN(f'eigV{number}', val=val)
        #
        NoClust = lattice.handle_no_clust(NoClust=args.NoClust)
        isingDictArgs = initialize_ising_dict_args(args, get_out_suffix(args), NoClust)
        #
        if NoClust is None:
            continue
        #
        isdy = IsingDynamics(lattice, **isingDictArgs)
        isdy.init_ising_dynamics(exName=isdy.run_id)
        match args.runlang:
            case 'C4'|'C1':
                lattice.export_ising_clust(exName=isdy.run_id,
                                        NoClust=NoClust)
                if args.runlang == 'C4':
                    lattice._export_eigV(number, exName=isdy.run_id, verbose=args.verbose)
        # Run in appropriate mode (equilibrium, SA, PT, or topo)
        rl = getattr(args, 'runlang', 'C1')
        if getattr(args, 'pt_mode', False):
            isdy.run(pt_mode=True, verbose=args.verbose)
        elif 'topo' in rl.lower():
            isdy.run(verbose=args.verbose)
        elif getattr(args, 'sa_mode', False):
            isdy.run(sa_mode=True, verbose=args.verbose)
        else:
            isdy.run(verbose=args.verbose)
        if getattr(args, 'save_results', False):
            _save_ising_results(isdy, lattice, args, _q)
        if args.remove_files:
            isdy.remove_run_c_files(remove_stderr=True)
            lattice.remove_exported_files()
