from typing import Any

from lrgsglib import ContactProcess, SignedGraph, join_non_empty

__all__ = [
    "initialize_contact_process_dict_args",
    "get_out_suffix",
    "run_contact_process",
    "clean_up_files",
]


def initialize_contact_process_dict_args(args: Any, out_suffix: str) -> dict[str, Any]:
    return dict(
        mu=args.mu,
        gamma=args.gamma,
        activation=args.activation,
        save_density=args.save_density,
        state_type=args.state_type,
        num_log_samples=args.num_log_samples,
        ic=args.init_cond,
        runlang=args.runlang,
        simpref=args.simpref,
        rndStr=args.rnd_str,
        out_suffix=out_suffix,
    )


def get_out_suffix(args: Any) -> str:
    return join_non_empty("_", args.init_cond, args.cell_type, args.out_suffix)


def run_contact_process(args: Any, signed_graph: SignedGraph) -> ContactProcess:
    contact_kwargs = initialize_contact_process_dict_args(args, get_out_suffix(args))
    cp = ContactProcess(signed_graph, **contact_kwargs)
    cp.init_contact_dynamics()
    cp.run(verbose=args.verbose, clean_export=args.remove_files)
    return cp


def clean_up_files(cp: ContactProcess, sg: SignedGraph, remove_stderr: bool = True) -> None:
    cp.remove_run_c_files(remove_stderr=remove_stderr)
    sg.remove_exported_files()
