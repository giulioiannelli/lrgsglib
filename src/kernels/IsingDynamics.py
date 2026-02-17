from lrgsglib import join_non_empty, IsingDynamics, SignedGraph, Any

__all__ = ["initialize_ising_dict_args", "get_out_suffix",
              "run_ising_dynamics", "clean_up_files"]

def initialize_ising_dict_args(args, out_suffix, NoClust):
    """Build IsingDynamics kwargs from parsed arguments."""
    base_args = dict(
        T=args.T,
        ic=args.init_cond,
        runlang=args.runlang,
        NoClust=NoClust,
        rndStr=args.randstr,
        freq=args.freq,
        out_suffix=out_suffix,
        thrmSTEP=args.thrmsteps
    )

    # Add SA parameters if sa_mode enabled
    if getattr(args, 'sa_mode', False):
        base_args.update(
            sa_enabled=True,
            T_init=args.T_init,
            T_final=args.T_final,
            cooling_schedule=args.cooling_schedule,
            cooling_rate=args.cooling_rate,
            steps_per_T=args.steps_per_T,
            n_temperatures=args.n_temperatures,
        )

    # Add PT parameters if pt_mode enabled
    if getattr(args, 'pt_mode', False):
        base_args.update(
            pt_enabled=True,
            n_replicas=args.n_replicas,
            T_min=args.T_min,
            T_max=args.T_max,
            T_ladder_type=args.T_ladder_type,
            steps_per_exchange=args.steps_per_exchange,
            n_exchanges=args.n_exchanges,
        )

    # Add topological algorithm parameters if topo runlang detected
    rl = getattr(args, 'runlang', 'C1')
    if 'topo' in rl.lower():
        base_args.update(
            topo_n_modes=getattr(args, 'topo_n_modes', 40),
            topo_sigma_init=getattr(args, 'topo_sigma_init', 0.15),
            topo_chunk_size=getattr(args, 'topo_chunk_size', 5000),
            topo_polish=getattr(args, 'topo_polish', True),
            topo_polish_sweeps=getattr(args, 'topo_polish_sweeps', 50),
            topo_tau=getattr(args, 'topo_tau', 1.0),
            topo_field_strength=getattr(args, 'topo_field_strength', 1.0),
        )
        # topo_fca delegates to SA, so ensure SA params are also set
        if 'topo_fca' in rl.lower() and not base_args.get('sa_enabled', False):
            base_args.update(
                sa_enabled=True,
                T_init=getattr(args, 'T_init', 10.0),
                T_final=getattr(args, 'T_final', 0.01),
                cooling_schedule=getattr(args, 'cooling_schedule', 'exponential'),
                cooling_rate=getattr(args, 'cooling_rate', 0.95),
                steps_per_T=getattr(args, 'steps_per_T', 100),
                n_temperatures=getattr(args, 'n_temperatures', 100),
            )

    return base_args

def get_out_suffix(args):
    return join_non_empty('_', args.init_cond, args.cell_type, args.out_suffix)

def run_ising_dynamics(
    args: Any,
    signed_graph: SignedGraph,
    out_suffix: str | None = None,
) -> Any:
    """
    Initialize and run the IsingDynamics simulation (equilibrium, SA, or PT mode).

    Parameters
    ----------
    args : Any
        Argument object containing Ising simulation parameters.
    signed_graph : SignedGraph
        Generic SignedGraph instance used for the dynamics.
    out_suffix : str or None
        Output suffix for file naming. If None, computed from args
        via :func:`get_out_suffix` (backward compat for L2D callers).

    Returns
    -------
    Any
        Executed IsingDynamics instance.
    """
    if out_suffix is None:
        out_suffix = get_out_suffix(args)
    isingDictArgs = initialize_ising_dict_args(args, out_suffix, 1)
    isdy = IsingDynamics(signed_graph, **isingDictArgs)
    isdy.init_ising_dynamics()

    rl = getattr(args, 'runlang', 'C1')
    if getattr(args, 'pt_mode', False):
        isdy.run(pt_mode=True)
    elif 'topo' in rl.lower():
        # Topological algorithms dispatch through the standard run() method
        # which detects topo_met/topo_fca from the runlang string
        isdy.run()
    elif getattr(args, 'sa_mode', False):
        isdy.run(sa_mode=True)
    else:
        isdy.run()

    return isdy


def clean_up_files(isdy: IsingDynamics, sg: SignedGraph) -> None:
    """
    Remove simulation-related temporary files.

    Parameters
    -----------
    isdy: IsingDynamics
        An IsingDynamics instance.
    sg: SignedGraph
        A SignedGraph instance.
    """
    if hasattr(isdy, 'remove_run_c_files'):
        isdy.remove_run_c_files()
    if hasattr(sg, 'remove_exported_files'):
        sg.remove_exported_files()