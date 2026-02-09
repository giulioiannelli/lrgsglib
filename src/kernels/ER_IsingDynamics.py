from typing import Any
from lrgsglib import IsingDynamics
from lrgsglib.graphs.nx import ErdosRenyi
from .ER import initialize_er_dict_args, prepare_erdos_renyi
from .IsingDynamics import (
    initialize_ising_dict_args, get_out_suffix,
    run_ising_dynamics, clean_up_files,
)
from parsers.shared import get_graph_engine

__all__ = ["run_simulation"]

MAX_ITER = 100_000


def run_simulation(
    args: Any,
    erDictArgs: dict[str, Any],
    isingDictArgs: dict[str, Any],
    number: int,
    remove_files: bool
) -> None:
    """
    Run Ising dynamics simulation on Erdos-Renyi graphs with retry logic.

    This function creates Erdos-Renyi graphs and runs Ising dynamics simulations,
    retrying if cluster formation fails. The retry mechanism is specific to ER graphs
    where not all realizations may form valid clusters.

    When ``--graph_engine gt`` is used, the eigenvector/cluster pipeline is
    skipped and the graph is created via ``prepare_erdos_renyi(args)``.

    Parameters
    ----------
    args : Any
        Argument object containing simulation configuration.
    erDictArgs : dict[str, Any]
        Keyword arguments for ErdosRenyi initialization.
    isingDictArgs : dict[str, Any]
        Keyword arguments for IsingDynamics initialization.
    number : int
        Eigenvector index to use for cluster initialization.
    remove_files : bool
        Whether to remove temporary files after simulation.
    """
    engine = get_graph_engine(args)

    if engine == "gt":
        _run_simulation_gt(args, remove_files)
    else:
        _run_simulation_nx(args, erDictArgs, isingDictArgs, number, remove_files)


def _run_simulation_gt(args: Any, remove_files: bool) -> None:
    """GT path: skip eigvec/cluster pipeline, use direct Ising dynamics."""
    for _ in range(args.number_of_averages):
        er = prepare_erdos_renyi(args)
        isdy = run_ising_dynamics(args, er)
        if remove_files:
            clean_up_files(isdy, er)


def _run_simulation_nx(
    args: Any,
    erDictArgs: dict[str, Any],
    isingDictArgs: dict[str, Any],
    number: int,
    remove_files: bool,
) -> None:
    """NX path: original eigenvector/cluster-based initialization with retry logic."""
    count = 0
    while count < args.number_of_averages and count < MAX_ITER:
        # Create new ER graph instance
        er = ErdosRenyi(**erDictArgs)
        er.flip_sel_edges(er.nwDict[args.cell_type]['G'])

        # Compute eigenvector and load on graph
        er.compute_k_eigvV(k=number + 1)
        er.load_eigV_on_graph(which=number, binarize=True)

        # Try to form clusters - may fail for some ER realizations
        try:
            er.make_clustersYN(f'eigV{number}', -1)
        except (ValueError, IndexError, KeyError) as exc:
            # Skip this realization if clustering fails
            continue

        # Run Ising dynamics on the successfully clustered graph
        isdy = IsingDynamics(er, **isingDictArgs)
        isdy.init_ising_dynamics()
        er.export_ising_clust()
        # run() with clean_export=True (default) handles file cleanup internally
        isdy.run(verbose=False, clean_export=remove_files)

        count += 1
