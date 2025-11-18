""" Contact process dynamics runner functions.
Command example to run EI dynamics on L2D contact process with batching of density files:
python lrgsglib/src/L2D_ContactProcess.py 64 0. -ac relu -na 200 -wd cptest -ga 0.493 -rl C1c -sp 1000 --randstr -sf 10
"""

from typing import Any
from pathlib import Path
import numpy as np
from lrgsglib.config.funcs import peq_fstr

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


_batch_densities = []
_last_saved_index = 0


def _process_EI_C1c(cp, args):
    """Batch process and aggregate density for EI+C1c, saving every (number_of_averages // save_frequency) runs."""
    global _batch_densities, _last_saved_index
    i = getattr(args, "_current_average", None)
    if i is None:
        raise ValueError("args._current_average must be set for aggregation filename.")
    na = getattr(args, "number_of_averages", None)
    sf = getattr(args, "save_frequency", getattr(args, "sf", None))
    if na is None or sf is None:
        sf = 1
        na = 1
    batch_size = max(1, na // sf)
    pdir = Path(getattr(cp.sg, "path_cntct", None) or getattr(cp.sg, "path_lrgsg", None) or cp.sg.path_data)
    dens_files = [
        f for f in pdir.glob("dens_*.bin") if cp.out_id in f.name or (getattr(cp, "rand_str", None) and cp.rand_str in f.name)
    ]
    if not dens_files:
        dens_files = sorted(
            [f for f in pdir.glob("dens_*.bin") if "_na=" not in f.name], key=lambda p: p.name
        )
    arr = None
    if dens_files:
        dens_file = dens_files[0]
        try:
            arr = np.fromfile(dens_file, dtype=np.float64)
        except Exception:
            arr = None
        try:
            dens_file.unlink()
        except Exception:
            pass
    if arr is not None:
        _batch_densities.append(arr)
    # Save at batch boundary or last run
    is_batch_end = (i % batch_size == 0) or (i == na)
    if is_batch_end and _batch_densities:
        # Build aggregate filename
        dynlabel = "gamma=" + str(getattr(args, "gamma", ""))
        out_suffix = get_out_suffix(args)
        agg_prefix = "_".join(filter(None, [peq_fstr(cp.sg.pflip), dynlabel, out_suffix]))
        sp_val = getattr(args, "sp", None)
        sp_token = f"_sp={int(sp_val)}" if sp_val is not None else ""
        agg_name = f"dens_{agg_prefix}{sp_token}_na={i}.bin"
        agg_path = pdir / agg_name
        concat = np.concatenate(_batch_densities)
        with open(agg_path, "wb") as f:
            concat.tofile(f)
        # Remove previous batch file
        if _last_saved_index > 0:
            prev_name = f"dens_{agg_prefix}{sp_token}_na={_last_saved_index}.bin"
            prev_path = pdir / prev_name
            if prev_path.exists():
                try:
                    prev_path.unlink()
                except Exception:
                    pass
        # Remove all earlier batch files
        for j in range(1, _last_saved_index):
            old_name = f"dens_{agg_prefix}{sp_token}_na={j}.bin"
            old_path = pdir / old_name
            if old_path.exists():
                try:
                    old_path.unlink()
                except Exception:
                    pass
        _batch_densities = []
        _last_saved_index = i
