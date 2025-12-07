from lrgsglib import Lattice2D

from .ContactProcessDynamics import run_contact_process, get_out_suffix, clean_up_files
from . import ContactProcessDynamics as cp_dyn
from .L2D import initialize_l2d_dict_args
from pathlib import Path
from dataclasses import dataclass
import numpy as np
from lrgsglib.config.funcs import peq_fstr


@dataclass
class RunInfo:
    out_id: str
    rand_str: str | None
    pflip: float
    path_cntct: str | None
    path_lrgsg: str | None
    path_data: str


def _post_process_EI_C1c(runinfos: list[RunInfo], args) -> None:
    """Aggregate and clean per-run data for EI+C1c runs (best-effort).

    This function concatenates per-run `dens_*.bin` files in the contact
    directory (order follows `runinfos`) into a single
    `dens_..._na=<n>.bin` and then removes per-run files that include the
    run's `out_id` or `rand_str`.
    """
    if not runinfos:
        return
    pdir = Path(runinfos[0].path_cntct or runinfos[0].path_lrgsg or runinfos[0].path_data)
    files = sorted([f for f in pdir.glob('dens_*.bin') if '_na=' not in f.name])
    rem = set(files)
    dens_arrays = []

    for ri in runinfos:
        found = next((f for f in rem if ri.out_id and ri.out_id in f.name), None)
        if not found and ri.rand_str:
            found = next((f for f in rem if ri.rand_str in f.name), None)
        if not found and rem:
            found = sorted(rem, key=lambda p: p.name)[0]
        if found:
            try:
                dens_arrays.append(np.fromfile(found, dtype=np.float64))
            except Exception:
                pass
            rem.discard(found)

        # remove any files in pdir that mention this run id / rand_str
        for f in list(pdir.glob('*')):
            if ri.out_id and ri.out_id in f.name:
                try:
                    f.unlink()
                except Exception:
                    pass
            elif ri.rand_str and ri.rand_str in f.name:
                try:
                    f.unlink()
                except Exception:
                    pass

    if dens_arrays:
        concat = np.concatenate(dens_arrays)
        dynlabel = f'gamma={getattr(args, "gamma", ""):.4g}'
        out_suffix = get_out_suffix(args)
        agg_prefix = '_'.join(filter(None, [peq_fstr(runinfos[0].pflip), dynlabel, out_suffix]))
        sp_val = getattr(args, 'sp', None)
        sp_token = f"_sp={int(sp_val)}" if sp_val is not None else ""
        out_path = pdir / f"dens_{agg_prefix}{sp_token}_na={int(args.number_of_averages)}.bin"
        concat.tofile(out_path)


def _post_process_default(cps: list, args) -> None:
    """Default post-processor for non-specialised backends.

    Currently this performs a best-effort per-run cleanup by delegating to
    `clean_up_files`. Keeping this as a separate function makes it easy to
    replace with aggregation logic for other dynamics/runlang combinations
    in the future.
    """
    for cp in cps:
        try:
            clean_up_files(cp, cp.sg, remove_stderr=True)
        except Exception:
            pass


def _prepare_lattice(args):
    lattice = Lattice2D(**initialize_l2d_dict_args(args))
    if lattice.init_nw_dict:
        lattice.flip_sel_edges(lattice.nwDict[args.cell_type]["G"])
    else:
        lattice.flip_random_fract_edges()
    return lattice


def _check_early_stopping(
    pdir: Path,
    agg_prefix: str,
    sp_token: str,
    last_saved_index: int,
    threshold: float,
    verbose: bool = False
) -> bool:
    """Check if early stopping condition is met.

    Args:
        pdir: Directory containing density files
        agg_prefix: Prefix for aggregated density filename
        sp_token: Sample preference token for filename
        last_saved_index: The last saved index (na value)
        threshold: Density threshold for early stopping
        verbose: Whether to print diagnostic messages

    Returns:
        True if early stopping condition is met, False otherwise
    """
    agg_name = f"dens_{agg_prefix}{sp_token}_na={last_saved_index}.bin"
    agg_path = pdir / agg_name

    if not agg_path.exists():
        return False

    try:
        dens_data = np.fromfile(agg_path, dtype=np.float64)
        if len(dens_data) > 0:
            # Reshape to (na, ns) where na is number of averages and ns is number of steps
            na = last_saved_index
            ns = len(dens_data) // na
            dens_reshaped = dens_data.reshape(na, ns)

            # Take last 10% of timesteps from each run
            M_per_run = max(1, int(0.1 * ns))
            last_M_per_run = dens_reshaped[:, -M_per_run:]
            avg_last_M = np.mean(last_M_per_run)

            if avg_last_M > threshold:
                if verbose:
                    print(f"""Early stopping triggered: last {M_per_run} 
                          timesteps per run (10% of {ns}) across {na} runs have 
                          average density {avg_last_M:.4f} > {threshold}""")
                return True
    except Exception as e:
        if verbose:
            print(f"Warning: Could not check early stopping condition: {e}")

    return False


def run_simulation(args):
    dynamics = args.dynamics.upper()
    runlang = args.runlang.upper()
    if dynamics == "EI" and runlang not in ("C1C", "C1D", "C1E", "C1F", "C1G", "PY"):
        raise NotImplementedError("Only C1c, C1d, C1e, C1f, C1g, and py backends are " \
        "currently wired for EI dynamics in L2D_ContactProcess.")
    if dynamics == "SIR":
        raise NotImplementedError("SIR-like contact-process wiring will " \
        "be implemented in a dedicated subprogram.")

    # Choose post-processing based on (dynamics, runlang).
    match (dynamics, runlang):
        case ("EI", "C1C") | ("EI", "C1D") | ("EI", "C1E") | ("EI", "C1F") | ("EI", "C1G"):
            from .ContactProcessDynamics import _process_EI_C1c, clean_up_files
            _process_cp = _process_EI_C1c
        case ("EI", "PY"):
            # Python backend uses default processor
            _process_cp = _post_process_default
        case _:
            raise NotImplementedError(f"""Post-processing for
                                      dynamics={dynamics} and runlang={runlang}
                                      is not implemented.""")

    # Detect existing aggregate to allow resuming.
    first_lattice = _prepare_lattice(args)
    pdir, agg_prefix, sp_token = cp_dyn._output_components(first_lattice, args)
    saved_idx = cp_dyn._latest_saved_index(pdir, agg_prefix, sp_token)
    cp_dyn._batch_densities = []
    cp_dyn._last_saved_index = saved_idx
    start_avg = saved_idx + 1
    if start_avg > args.number_of_averages:
        if args.verbose:
            print(f"Density file already has na={saved_idx}; nothing to do.")
        return

    # Early stopping parameters
    early_stop_threshold = args.early_stop_density_threshold
    min_runs_before_check = 20

    for i in range(start_avg, args.number_of_averages + 1):
        lattice = first_lattice if i == start_avg else _prepare_lattice(args)
        cp = run_contact_process(args, lattice)
        setattr(args, '_current_average', i)
        _process_cp(cp, args)
        try:
            clean_up_files(cp, cp.sg, remove_stderr=True)
        except Exception:
            pass

        # Check early stopping condition after minimum runs
        if early_stop_threshold is not None and i >= min_runs_before_check:
            if _check_early_stopping(pdir, agg_prefix, sp_token, 
                                     cp_dyn._last_saved_index, 
                                     early_stop_threshold, args.verbose):
                if args.verbose:
                    print(f"Stopping at run {i}")
                break
