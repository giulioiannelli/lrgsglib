from typing import Any

from lrgsglib import ContactProcessEI, ContactProcessSIR, SignedGraph, join_non_empty

__all__ = [
    "initialize_contact_process_dict_args",
    "get_out_suffix",
    "run_contact_process",
    "clean_up_files",
]


def _normalize_dynamics(args: Any) -> str:
    if hasattr(args, "dynamics") and args.dynamics is not None:
        return str(args.dynamics).upper()
    runlang = str(getattr(args, "runlang", "")).upper()
    return "EI" if runlang.startswith("C1") else "SIR"


def _validate_backend_choice(dynamics: str, runlang: str) -> None:
    if dynamics == "EI" and not runlang.startswith("C1"):
        raise ValueError("Excitation/inhibition dynamics require a C1* runlang (e.g., C1c).")
    if dynamics == "SIR" and runlang.startswith("C1"):
        raise ValueError("SIR dynamics do not map to the C1* contact-process kernels.")


def initialize_contact_process_dict_args(args: Any, out_suffix: str) -> dict[str, Any]:
    dynamics = _normalize_dynamics(args)
    runlang = str(args.runlang).upper()
    _validate_backend_choice(dynamics, runlang)

    common: dict[str, Any] = dict(
        save_density=args.save_density,
        state_type=args.state_type,
        ic=args.init_cond,
        runlang=args.runlang,
        simpref=args.simpref,
        rndStr=args.randstr,
        out_suffix=out_suffix,
    )

    if dynamics == "EI":
        if args.gamma is None:
            raise ValueError("Argument gamma is required for excitation/inhibition dynamics.")
        common.update(
            gamma=args.gamma,
            activation=args.activation,
            num_log_samples=args.num_log_samples,
        )
    else:
        if args.mu is None:
            raise ValueError("Argument mu is required for SIR-style dynamics.")
        common.update(mu=args.mu)
    return common


def get_out_suffix(args: Any) -> str:
    return join_non_empty("_", args.init_cond, args.cell_type, args.out_suffix)


def run_contact_process(args: Any, signed_graph: SignedGraph):
    contact_kwargs = initialize_contact_process_dict_args(args, get_out_suffix(args))
    dynamics = _normalize_dynamics(args)
    cp_cls = ContactProcessEI if dynamics == "EI" else ContactProcessSIR
    cp = cp_cls(signed_graph, **contact_kwargs)
    cp.init_contact_dynamics()
    cp.run(verbose=args.verbose, clean_export=args.remove_files)
    return cp


def clean_up_files(cp: Any, sg: SignedGraph, remove_stderr: bool = True) -> None:
    cp.remove_run_c_files(remove_stderr=remove_stderr)
    sg.remove_exported_files()
