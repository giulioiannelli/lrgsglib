""" Contact process dynamics runner functions.
Command example to run EI dynamics on L2D contact process with batching of density files:
python lrgsglib/src/L2D_ContactProcess.py 64 0. -ac relu -na 200 -wd cptest -ga 0.493 -rl C1c -sp 1000 --randstr -sf 10

Supported C1 backends for EI dynamics: C1c, C1d, C1e, C1f, C1g
- C1c: Log-spaced density sampling with standard MC sweeps
- C1d: Adaptive frontier optimization for low-density regimes
- C1e: Cached-lambda updates without frontier (log-spaced sampling)
- C1f: Gillespie-style event-driven loop over frontier
- C1g: C1e with configuration snapshots at log-spaced intervals
"""

from typing import Any
from pathlib import Path
import re
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


def _output_components(sg: Any, args: Any) -> tuple[Path, str, str]:
    """Return (pdir, agg_prefix, sp_token) for density files."""
    pdir = Path(getattr(sg, "path_cntct", None) or getattr(sg, "path_lrgsg", None) or sg.path_data)
    dynlabel = f'gamma={getattr(args, "gamma", ""):.3g}'
    out_suffix = get_out_suffix(args)
    agg_prefix = "_".join(filter(None, [peq_fstr(sg.pflip), dynlabel, out_suffix]))
    sp_val = getattr(args, "sp", None)
    sp_token = f"_sp={int(sp_val)}" if sp_val is not None else ""
    return pdir, agg_prefix, sp_token


def _latest_saved_index(pdir: Path, agg_prefix: str, sp_token: str) -> int:
    """Find the largest `_na=<n>` aggregate file already on disk."""
    pattern = re.compile(rf"^dens_{re.escape(agg_prefix)}{re.escape(sp_token)}_na=(\d+)\.bin$")
    latest = 0
    for f in pdir.glob(f"dens_{agg_prefix}{sp_token}_na=*.bin"):
        m = pattern.match(f.name)
        if m:
            try:
                latest = max(latest, int(m.group(1)))
            except ValueError:
                pass
    return latest


def _process_EI_C1c(cp, args):
    """Batch process and aggregate density for EI+C1c.

    Saves a cumulative file at every batch boundary with the filename
    suffix `_na=<current>` so early stops can be resumed. Only the latest
    `_na=<current>` file is kept.
    """
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
    pdir, agg_prefix, sp_token = _output_components(cp.sg, args)
    if _last_saved_index == 0:
        _last_saved_index = _latest_saved_index(pdir, agg_prefix, sp_token)
    # Collect per-run density for this batch.
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
        # Load previously aggregated data if present (for resume or prior batch).
        prev_data = None
        if _last_saved_index > 0:
            prev_name = f"dens_{agg_prefix}{sp_token}_na={_last_saved_index}.bin"
            prev_path = pdir / prev_name
            if prev_path.exists():
                try:
                    prev_data = np.fromfile(prev_path, dtype=np.float64)
                except Exception:
                    prev_data = None
        if prev_data is None:
            prev_data = np.array([], dtype=np.float64)

        concat = np.concatenate(_batch_densities)
        new_data = np.concatenate([prev_data, concat])

        agg_name = f"dens_{agg_prefix}{sp_token}_na={i}.bin"
        agg_path = pdir / agg_name
        with open(agg_path, "wb") as f:
            new_data.tofile(f)

        # Remove older aggregate files to keep only the latest progress marker.
        for f in pdir.glob(f"dens_{agg_prefix}{sp_token}_na=*.bin"):
            m = re.search(r"_na=(\d+)\.bin$", f.name)
            if m:
                idx = int(m.group(1))
                if idx < i:
                    try:
                        f.unlink()
                    except Exception:
                        pass

        _batch_densities = []
        _last_saved_index = i
