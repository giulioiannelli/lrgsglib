from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from lrgsglib import join_non_empty, IsingDynamics, SignedGraph
from .generic import find_existing_averages

__all__ = [
    "initialize_ising_dict_args",
    "get_out_suffix",
    "run_ising_dynamics",
    "clean_up_files",
    # Checkpoint utilities
    "build_ising_stem",
    "load_ising_checkpoint",
    "save_ising_checkpoint",
    "find_completed_quenches",
    "extract_ene_data",
    "extract_magn_data",
    "extract_sout_data",
    "extract_coeffs_data",
    "new_accumulator",
    "rebuild_accumulator",
    "append_run",
    "finalize_accumulator",
    "should_checkpoint",
    "resolve_save_flags",
]

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
        thrmSTEP=args.thrmsteps,
        eqSTEP=getattr(args, 'eqstep', 20),
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
        # topo_fca quench uses n_temperatures*steps_per_T for total sweeps
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


# ── Checkpoint utilities ──────────────────────────────────────────────

# Keys whose arrays are stacked along axis-0 (one row per thermal run).
_STACKABLE_ENE = ("ene", "sa_energy")
_STACKABLE_MAGN = ("magn", "sa_magn")
_STACKABLE_SOUT = ("s_t",)
_STACKABLE_COEFFS = (
    "topo_met_coeffs", "topo_met_best_spins", "topo_met_best_energy",
)
# Keys stored once (shared across thermal runs).
_SHARED_KEYS = ("sa_temps",)
# All stackable keys (union).
_ALL_STACKABLE = (
    _STACKABLE_ENE + _STACKABLE_MAGN + _STACKABLE_SOUT + _STACKABLE_COEFFS
)


def resolve_save_flags(args: Any) -> dict[str, bool]:
    """Resolve which observables to save from CLI flags.

    ``--save-all`` overrides individual flags to True.

    Returns
    -------
    dict[str, bool]
        Keys: ``"ene"``, ``"magn"``, ``"sout"``, ``"coeffs"``.
    """
    do_save = getattr(args, "save_results", False)
    if not do_save:
        return {"ene": False, "magn": False, "sout": False, "coeffs": False}
    save_all = getattr(args, "save_all", False)
    rl = getattr(args, "runlang", "").lower()
    return {
        "ene": save_all or getattr(args, "save_ene", True),
        "magn": save_all or getattr(args, "save_magn", True),
        "sout": save_all or getattr(args, "save_sout", False),
        "coeffs": (save_all or True) if "topo_met" in rl else False,
    }


def build_ising_stem(
    prefix: str,
    rl: str,
    lattice: Any,
    quench_idx: int,
    *,
    is_l3d: bool = False,
) -> tuple[str, Path]:
    """Build filename stem and save directory for an observable file.

    Parameters
    ----------
    prefix : str
        Observable prefix (``"ene"``, ``"magn"``, ``"sout"``, ``"coeffs"``).
    rl : str
        Runlang string (e.g. ``"pb_sa"``).
    lattice : SignedGraph
        Lattice instance (used for side, pflip, geo, path_sgdata).
    quench_idx : int
        Quench realization index.
    is_l3d : bool
        If True extract side from ``lattice.dim[0]`` (3D), else
        ``lattice.side1`` (2D).

    Returns
    -------
    tuple[str, Path]
        ``(stem, save_dir)`` where stem is e.g.
        ``"ene_pb_sa_L=8_sc_p=0.3_q=000"`` and save_dir is the target
        directory.
    """
    if is_l3d:
        dim = getattr(lattice, "dim", (0,))
        side = dim[0] if isinstance(dim, tuple) else dim
    else:
        side = getattr(lattice, "side1", 0)
    pflip = getattr(lattice, "pflip", 0.0)
    geo = getattr(lattice, "geo", "sqr" if not is_l3d else "sc")

    is_topo = "topo" in rl.lower()
    n_modes = getattr(lattice, "_topo_n_modes", 0)
    # topo_n_modes may be on the args, not lattice; caller should set
    # lattice._topo_n_modes before calling if needed.  Alternatively
    # accept it from the IsingDynamics instance.  We fall back to 0.
    mode_tag = f"_m={n_modes}" if is_topo and n_modes else ""

    stem = (
        f"{prefix}_{rl}_L={side}_{geo}_p={pflip:.3g}"
        f"{mode_tag}_q={quench_idx:03d}"
    )
    save_dir = Path(lattice.path_sgdata) / "ising_results"
    return stem, save_dir


def load_ising_checkpoint(
    save_dir: Path,
    stem: str,
) -> tuple[dict[str, Any] | None, int, Path | None]:
    """Load the most recent checkpoint matching *stem*.

    Searches for files ``{stem}_nt=*.npz`` and returns the one with the
    highest ``nt`` value.

    Returns
    -------
    tuple[dict | None, int, Path | None]
        ``(data, n_thermal_done, checkpoint_path)``.
        ``(None, 0, None)`` when no checkpoint exists.
    """
    if not save_dir.exists():
        return None, 0, None

    pattern = f"{stem}_nt=*.npz"
    n_done = find_existing_averages(
        save_dir, pattern, count_token="nt", find_max=True,
    )
    if n_done == 0:
        return None, 0, None

    # Find the actual file with nt=n_done
    for fpath in save_dir.glob(pattern):
        parts = [p for p in fpath.stem.split("_") if p.startswith("nt=")]
        if parts and int(parts[0].split("=")[1]) == n_done:
            data = dict(np.load(fpath, allow_pickle=True))
            return data, n_done, fpath

    return None, 0, None


def find_completed_quenches(
    save_dir: Path,
    stem_no_q: str,
    target_nt: int,
    *,
    prefixes: list[str] | None = None,
) -> set[int]:
    """Return quench indices that are fully complete.

    A quench is complete when **all** requested observable prefixes
    have ``nt == target_nt``.

    Parameters
    ----------
    save_dir : Path
        Directory containing checkpoint files.
    stem_no_q : str
        The filename stem with ``_q=*`` as a wildcard placeholder.
        E.g. ``"ene_pb_sa_L=8_sc_p=0.3"``.
    target_nt : int
        Required number of thermal runs for completion.
    prefixes : list[str] or None
        Observable prefixes to check (e.g. ``["ene", "magn"]``).
        If None, checks only for the prefix already in *stem_no_q*.

    Returns
    -------
    set[int]
        Quench indices where all observables have ``nt == target_nt``.
    """
    if not save_dir.exists():
        return set()

    # Discover all q indices from the first prefix
    pattern = f"{stem_no_q}_q=*_nt={target_nt}.npz"
    q_sets: list[set[int]] = []

    if prefixes is None:
        prefixes_to_check = [stem_no_q.split("_")[0]]
    else:
        prefixes_to_check = prefixes

    for pfx in prefixes_to_check:
        # Replace the leading prefix in stem_no_q
        parts = stem_no_q.split("_", 1)
        if len(parts) > 1:
            stem_pfx = f"{pfx}_{parts[1]}"
        else:
            stem_pfx = pfx
        pat = f"{stem_pfx}_q=*_nt={target_nt}.npz"
        qs: set[int] = set()
        for fpath in save_dir.glob(pat):
            q_parts = [
                p for p in fpath.stem.split("_") if p.startswith("q=")
            ]
            if q_parts:
                qs.add(int(q_parts[0].split("=")[1]))
        q_sets.append(qs)

    if not q_sets:
        return set()
    # Intersection: complete only if all prefixes are done
    return set.intersection(*q_sets)


# ── Data extraction (one run → dict) ─────────────────────────────────

def extract_ene_data(isdy: Any) -> dict[str, Any]:
    """Extract energy arrays from a single IsingDynamics run.

    For SA/TFCA: ``sa_energy`` + ``sa_temps`` (shared).
    For equilibrium/topo_met: ``ene``.
    """
    d: dict[str, Any] = {}
    if hasattr(isdy, "sa_energy") and isdy.sa_energy is not None:
        d["sa_energy"] = np.asarray(isdy.sa_energy)
        if hasattr(isdy, "sa_temps") and isdy.sa_temps is not None:
            d["sa_temps"] = np.asarray(isdy.sa_temps)
    if hasattr(isdy, "ene") and isdy.ene is not None:
        d["ene"] = np.asarray(isdy.ene)
    return d


def extract_magn_data(isdy: Any) -> dict[str, Any]:
    """Extract magnetization arrays from a single run."""
    d: dict[str, Any] = {}
    if hasattr(isdy, "sa_magn") and isdy.sa_magn is not None:
        d["sa_magn"] = np.asarray(isdy.sa_magn)
    if hasattr(isdy, "magn") and isdy.magn is not None:
        d["magn"] = np.asarray(isdy.magn)
    return d


def extract_sout_data(isdy: Any) -> dict[str, Any]:
    """Extract spin snapshot arrays from a single run."""
    d: dict[str, Any] = {}
    if hasattr(isdy, "s_t") and isdy.s_t is not None and len(isdy.s_t) > 0:
        # s_t is a list of 1D arrays → stack to 2D (n_snapshots, N)
        d["s_t"] = np.stack(isdy.s_t, axis=0)
    return d


def extract_coeffs_data(isdy: Any) -> dict[str, Any]:
    """Extract topo_met coefficients and best-spin data."""
    d: dict[str, Any] = {}
    if hasattr(isdy, "topo_met_coeffs") and isdy.topo_met_coeffs is not None:
        d["topo_met_coeffs"] = np.asarray(isdy.topo_met_coeffs)
    if hasattr(isdy, "topo_met_best_energy"):
        d["topo_met_best_energy"] = np.float64(isdy.topo_met_best_energy)
    if hasattr(isdy, "topo_met_best_spins") and isdy.topo_met_best_spins is not None:
        d["topo_met_best_spins"] = np.asarray(isdy.topo_met_best_spins)
    return d


# ── Accumulator management ───────────────────────────────────────────

def _build_metadata(
    lattice: Any,
    args: Any,
    quench_idx: int,
    *,
    is_l3d: bool = False,
) -> dict[str, Any]:
    """Build the scalar metadata dict for an NPZ file."""
    rl = getattr(args, "runlang", "C1")
    if is_l3d:
        dim = getattr(lattice, "dim", (0,))
        side = dim[0] if isinstance(dim, tuple) else dim
    else:
        side = getattr(lattice, "side1", 0)
    pflip = getattr(lattice, "pflip", 0.0)
    geo = getattr(lattice, "geo", "sqr" if not is_l3d else "sc")
    G = lattice.gr["G"]
    n_edges = G.number_of_edges()

    meta: dict[str, Any] = {
        "runlang": rl,
        "side": side,
        "pflip": pflip,
        "geo": geo,
        "N": G.number_of_nodes(),
        "n_edges": n_edges,
        "quench_idx": quench_idx,
    }
    if is_l3d:
        meta["dim"] = dim

    # Topo params
    if "topo" in rl.lower():
        for key in (
            "topo_n_modes", "topo_sigma_init", "topo_chunk_size",
            "topo_polish", "topo_polish_sweeps", "topo_tau",
            "topo_field_strength",
        ):
            meta[key] = getattr(args, key, 0)
    return meta


def new_accumulator(
    lattice: Any,
    args: Any,
    quench_idx: int,
    *,
    is_l3d: bool = False,
) -> dict[str, Any]:
    """Create a fresh empty accumulator.

    Internal layout::

        {"_meta": {...}, "_lists": {}, "_shared": {}}

    ``_lists`` values are ``list[np.ndarray]`` (one entry per thermal
    run).  ``_shared`` stores arrays that are identical across runs.
    """
    return {
        "_meta": _build_metadata(
            lattice, args, quench_idx, is_l3d=is_l3d,
        ),
        "_lists": {},
        "_shared": {},
    }


def rebuild_accumulator(data: dict[str, Any]) -> dict[str, Any]:
    """Reconstruct an accumulator from a loaded NPZ checkpoint.

    Splits 2-D stacked arrays back into lists of 1-D (or lower-dim)
    arrays so that :func:`append_run` can continue adding rows.
    """
    n_t = int(data.get("n_thermal", 1))
    meta: dict[str, Any] = {}
    lists: dict[str, list[np.ndarray]] = {}
    shared: dict[str, np.ndarray] = {}

    for key, val in data.items():
        arr = np.asarray(val)
        if key == "n_thermal":
            continue
        if key in _ALL_STACKABLE and arr.ndim >= 1 and arr.shape[0] == n_t:
            # Stacked array → split back to list
            lists[key] = [arr[i] for i in range(n_t)]
        elif key in _SHARED_KEYS:
            shared[key] = arr
        else:
            meta[key] = val

    return {"_meta": meta, "_lists": lists, "_shared": shared}


def append_run(accumulated: dict[str, Any], single: dict[str, Any]) -> None:
    """Append one thermal run's data to the accumulator (in-place)."""
    lists = accumulated["_lists"]
    shared = accumulated["_shared"]
    for key, val in single.items():
        if key in _SHARED_KEYS:
            if key not in shared:
                shared[key] = np.asarray(val)
        elif key in _ALL_STACKABLE:
            if key not in lists:
                lists[key] = []
            lists[key].append(np.asarray(val))
        # scalar topo metadata is already in _meta


def finalize_accumulator(
    accumulated: dict[str, Any],
    n_thermal: int,
) -> dict[str, Any]:
    """Convert accumulator to a flat dict suitable for ``np.savez``."""
    result = dict(accumulated["_meta"])
    result["n_thermal"] = n_thermal
    for key, arr_list in accumulated["_lists"].items():
        if arr_list:
            result[key] = np.stack(arr_list, axis=0)
    for key, arr in accumulated["_shared"].items():
        result[key] = arr
    return result


# ── Checkpoint I/O ────────────────────────────────────────────────────

def save_ising_checkpoint(
    data: dict[str, Any],
    save_dir: Path,
    stem: str,
    n_thermal: int,
    *,
    old_path: Path | None = None,
) -> Path:
    """Write accumulated data to an NPZ checkpoint.

    Removes the previous checkpoint if its path differs from the new
    one (atomic rename pattern).

    Returns
    -------
    Path
        Path to the newly saved checkpoint.
    """
    save_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{stem}_nt={n_thermal}.npz"
    new_path = save_dir / fname
    np.savez_compressed(new_path, **data)
    if old_path and old_path != new_path and old_path.exists():
        old_path.unlink()
    return new_path


def should_checkpoint(
    n_done: int, freq: int, target: int,
) -> bool:
    """Return True if a checkpoint should be written now.

    Parameters
    ----------
    n_done : int
        Number of thermal runs completed so far.
    freq : int
        Save frequency (0 means save only at the end).
    target : int
        Total number of thermal runs requested.
    """
    if n_done == target:
        return True
    if freq <= 0:
        return False
    return n_done % freq == 0