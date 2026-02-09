from typing import Any
from lrgsglib import IsingDynamics
from lrgsglib.graphs.nx import ErdosRenyi
from .ER import initialize_er_dict_args
from .IsingDynamics import initialize_ising_dict_args, get_out_suffix

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

    Notes
    -----
    The function will retry up to MAX_ITER times if cluster formation fails
    (make_clustersYN raises an exception). This is ER-specific behavior as
    not all ER realizations may have sufficient structure for clustering.

    Examples
    --------
    >>> class Args:
    ...     number_of_averages = 10
    ...     cell_type = 'rand'
    ...     thrmsteps = 20
    >>> erDictArgs = {'n': 100, 'p': 0.1, 'sgpathn': '', 'pflip': 0.2, 'init_nw_dict': True}
    >>> isingDictArgs = {'T': 1.0, 'ic': 'ground_state_0', 'runlang': 'C1', 'NoClust': 1}
    >>> run_simulation(Args(), erDictArgs, isingDictArgs, 0, True)  # doctest: +SKIP
    """
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
