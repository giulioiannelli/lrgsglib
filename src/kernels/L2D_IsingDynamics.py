"""Kernel for Ising dynamics on 2D lattices.

Supports three execution paths:
- GT engine (skip eigvec/cluster pipeline)
- NX direct (no ground-state init)
- NX ground-state (eigenvector/cluster-based init)

All paths share the same checkpoint/accumulation logic via
:mod:`kernels.IsingDynamics`.
"""

from __future__ import annotations

from typing import Any

from lrgsglib import *
from .L2D import *
from .IsingDynamics import (
    initialize_ising_dict_args,
    get_out_suffix,
    run_ising_dynamics,
    clean_up_files,
    resolve_save_flags,
    build_ising_stem,
    load_ising_checkpoint,
    save_ising_checkpoint,
    find_completed_quenches,
    extract_ene_data,
    extract_magn_data,
    extract_sout_data,
    extract_coeffs_data,
    new_accumulator,
    rebuild_accumulator,
    append_run,
    finalize_accumulator,
    should_checkpoint,
)
from parsers.shared import get_graph_engine

# Map observable prefix → extractor function
_EXTRACTORS = {
    "ene": extract_ene_data,
    "magn": extract_magn_data,
    "sout": extract_sout_data,
    "coeffs": extract_coeffs_data,
}


# ── Checkpoint helpers (shared by all L2D paths) ─────────────────────

def _setup_obs_state(
    save_flags: dict[str, bool],
    rl: str,
    lattice: Any,
    quench_idx: int,
    args: Any,
) -> dict[str, dict[str, Any]]:
    """Build per-observable accumulator state for one quench.

    Returns a dict keyed by prefix ("ene", "magn", ...) with entries:
    ``stem``, ``dir``, ``acc``, ``n_done``, ``old_path``.
    """
    if "topo" in rl.lower():
        lattice._topo_n_modes = getattr(args, "topo_n_modes", 0)

    obs_state: dict[str, dict[str, Any]] = {}
    for prefix, enabled in save_flags.items():
        if not enabled:
            continue
        stem, save_dir = build_ising_stem(
            prefix, rl, lattice, quench_idx, is_l3d=False,
        )
        prior, n_done, old_path = load_ising_checkpoint(save_dir, stem)
        acc = (
            rebuild_accumulator(prior)
            if prior is not None
            else new_accumulator(lattice, args, quench_idx, is_l3d=False)
        )
        obs_state[prefix] = {
            "stem": stem,
            "dir": save_dir,
            "acc": acc,
            "n_done": n_done,
            "old_path": old_path,
        }
    return obs_state


def _do_checkpoint(
    obs_state: dict[str, dict[str, Any]],
    current_nt: int,
    save_freq: int,
    n_thermal: int,
) -> None:
    """Write checkpoints for all observables if it's time."""
    if not should_checkpoint(current_nt, save_freq, n_thermal):
        return
    for state in obs_state.values():
        final = finalize_accumulator(state["acc"], current_nt)
        new_path = save_ising_checkpoint(
            final, state["dir"], state["stem"],
            current_nt, old_path=state["old_path"],
        )
        state["old_path"] = new_path


def _accumulate_run(
    obs_state: dict[str, dict[str, Any]],
    isdy: Any,
) -> None:
    """Extract and append one run's data to all active accumulators."""
    for prefix, state in obs_state.items():
        single = _EXTRACTORS[prefix](isdy)
        append_run(state["acc"], single)


# ── Main entry point ─────────────────────────────────────────────────

def run_simulation(args):
    engine = get_graph_engine(args)

    if engine == "gt":
        _run_simulation_gt(args)
    else:
        _run_simulation_nx(args)


def _run_simulation_gt(args):
    """GT path: skip eigvec/cluster pipeline, use direct Ising dynamics."""
    n_thermal = getattr(args, "n_thermal", 1)
    save_freq = getattr(args, "save_frequency", 0)
    save_flags = resolve_save_flags(args)
    any_save = any(save_flags.values())

    for _q in range(1, args.number_of_averages + 1):
        lattice = prepare_lattice(args)

        obs_state: dict[str, dict[str, Any]] = {}
        n_done = 0
        if any_save:
            rl = getattr(args, "runlang", "C1")
            obs_state = _setup_obs_state(
                save_flags, rl, lattice, _q, args,
            )
            n_done = (
                min(s["n_done"] for s in obs_state.values())
                if obs_state
                else 0
            )

        for _t in range(n_done, n_thermal):
            isdy = run_ising_dynamics(args, lattice)
            if obs_state:
                _accumulate_run(obs_state, isdy)
                _do_checkpoint(obs_state, _t + 1, save_freq, n_thermal)
            if args.remove_files:
                clean_up_files(isdy, lattice)


def _run_simulation_nx(args):
    """NX path: eigenvector/cluster-based init for ground_state, direct for others."""
    n_thermal = getattr(args, "n_thermal", 1)
    save_freq = getattr(args, "save_frequency", 0)
    save_flags = resolve_save_flags(args)
    any_save = any(save_flags.values())
    ic_gs = (
        args.init_cond.startswith("ground_state")
        or args.init_cond.startswith("gs")
    )

    if not ic_gs:
        _run_nx_direct(args, n_thermal, save_freq, save_flags, any_save)
    else:
        _run_nx_ground_state(
            args, n_thermal, save_freq, save_flags, any_save,
        )


def _run_nx_direct(args, n_thermal, save_freq, save_flags, any_save):
    """NX direct path: no eigenvector/cluster pipeline needed."""
    for _q in range(1, args.number_of_averages + 1):
        lattice = Lattice2D(**initialize_l2d_dict_args(args))
        if lattice.init_nw_dict:
            lattice.flip_sel_edges(lattice.nwDict[args.cell_type]["G"])
        else:
            lattice.flip_random_fract_edges()

        obs_state: dict[str, dict[str, Any]] = {}
        n_done = 0
        if any_save:
            rl = getattr(args, "runlang", "C1")
            obs_state = _setup_obs_state(
                save_flags, rl, lattice, _q, args,
            )
            n_done = (
                min(s["n_done"] for s in obs_state.values())
                if obs_state
                else 0
            )

        for _t in range(n_done, n_thermal):
            isdy = run_ising_dynamics(args, lattice)
            if obs_state:
                _accumulate_run(obs_state, isdy)
                _do_checkpoint(obs_state, _t + 1, save_freq, n_thermal)
            if args.remove_files:
                clean_up_files(isdy, lattice)


def _run_nx_ground_state(
    args, n_thermal, save_freq, save_flags, any_save,
):
    """NX ground-state path: eigenvector/cluster-based initialization."""
    number = int(args.init_cond.split("_")[-1])
    val = ConditionalPartitioning(args.val)

    for _q in range(1, args.number_of_averages + 1):
        lattice = Lattice2D(**initialize_l2d_dict_args(args))
        if lattice.init_nw_dict:
            lattice.flip_sel_edges(lattice.nwDict[args.cell_type]["G"])
        else:
            lattice.flip_random_fract_edges()
        lattice.compute_k_eigvV(number + 1)
        lattice.load_eigV_on_graph(number, binarize=True)
        lattice.make_clustersYN(f"eigV{number}", val=val)

        NoClust = lattice.handle_no_clust(NoClust=args.NoClust)
        if NoClust is None:
            continue

        isingDictArgs = initialize_ising_dict_args(
            args, get_out_suffix(args), NoClust,
        )
        isdy = IsingDynamics(lattice, **isingDictArgs)
        isdy.init_ising_dynamics(exName=isdy.run_id)
        match args.runlang:
            case "C4" | "C1":
                lattice.export_ising_clust(
                    exName=isdy.run_id, NoClust=NoClust,
                )
                if args.runlang == "C4":
                    lattice._export_eigV(
                        number, exName=isdy.run_id,
                        verbose=args.verbose,
                    )

        obs_state: dict[str, dict[str, Any]] = {}
        n_done = 0
        if any_save:
            rl = getattr(args, "runlang", "C1")
            obs_state = _setup_obs_state(
                save_flags, rl, lattice, _q, args,
            )
            n_done = (
                min(s["n_done"] for s in obs_state.values())
                if obs_state
                else 0
            )

        for _t in range(n_done, n_thermal):
            # Run in appropriate mode
            rl = getattr(args, "runlang", "C1")
            if getattr(args, "pt_mode", False):
                isdy.run(pt_mode=True, verbose=args.verbose)
            elif "topo" in rl.lower():
                isdy.run(verbose=args.verbose)
            elif getattr(args, "sa_mode", False):
                isdy.run(sa_mode=True, verbose=args.verbose)
            else:
                isdy.run(verbose=args.verbose)

            if obs_state:
                _accumulate_run(obs_state, isdy)
                _do_checkpoint(obs_state, _t + 1, save_freq, n_thermal)
            if args.remove_files:
                isdy.remove_run_c_files(remove_stderr=True)
                lattice.remove_exported_files()
