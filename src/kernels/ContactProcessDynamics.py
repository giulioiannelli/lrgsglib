from typing import Any

from lrgsglib import ContactProcessEI, ContactProcessSIR, SignedGraph, join_non_empty

__all__ = [
    "initialize_contact_process_dict_args",
    "get_out_suffix",
    "run_contact_process",
    "clean_up_files",
]


def initialize_contact_process_dict_args(args: Any, out_suffix: str) -> dict[str, Any]:
    common: dict[str, Any] = dict(
        save_density=args.save_density,
        state_type=args.state_type,
        ic=args.init_cond,
        runlang=args.runlang,
        simpref=args.simpref,
        rndStr=args.rnd_str,
        out_suffix=out_suffix,
    )
    if str(args.runlang).upper().startswith("C1"):
        common.update(
            gamma=args.gamma,
            activation=args.activation,
            num_log_samples=args.num_log_samples,
        )
    else:
        common.update(mu=args.mu)
    return common


def get_out_suffix(args: Any) -> str:
    return join_non_empty("_", args.init_cond, args.cell_type, args.out_suffix)


def run_contact_process(args: Any, signed_graph: SignedGraph):
    contact_kwargs = initialize_contact_process_dict_args(args, get_out_suffix(args))
    cp_cls = ContactProcessEI if str(args.runlang).upper().startswith("C1") else ContactProcessSIR
    cp = cp_cls(signed_graph, **contact_kwargs)
    cp.init_contact_dynamics()
    cp.run(verbose=args.verbose, clean_export=args.remove_files)
    return cp


def clean_up_files(cp: Any, sg: SignedGraph, remove_stderr: bool = True) -> None:
    cp.remove_run_c_files(remove_stderr=remove_stderr)
    sg.remove_exported_files()
