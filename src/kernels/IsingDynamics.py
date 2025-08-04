from lrgsglib import join_non_empty, IsingDynamics, SignedGraph, Any

__all__ = ["initialize_ising_dict_args", "get_out_suffix",
              "run_ising_dynamics", "clean_up_files"]

def initialize_ising_dict_args(args, out_suffix, NoClust):
    return dict(
        T=args.T, 
        ic=args.init_cond, 
        runlang=args.runlang, 
        NoClust=NoClust, 
        rndStr=args.rnd_str,
        freq=args.freq,
        out_suffix=out_suffix,
        thrmSTEP=args.thrmsteps
    )

def get_out_suffix(args):
    return join_non_empty('_', args.init_cond, args.cell_type, args.out_suffix)

def run_ising_dynamics(args: Any, signed_graph: SignedGraph) -> Any:
    """
    Initialize and run the IsingDynamics simulation.

    Parameters
    ----------
    args : Any
        Argument object containing Ising simulation parameters.
    signed_graph : SignedGraph
        Generic SignedGraph instance used for the dynamics.

    Returns
    -------
    Any
        Executed IsingDynamics instance.
    """
    isingDictArgs = initialize_ising_dict_args(args, get_out_suffix(args), 1)
    isdy = IsingDynamics(signed_graph, **isingDictArgs)
    isdy.init_ising_dynamics()
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
    isdy.remove_run_c_files()
    sg.remove_exported_files()