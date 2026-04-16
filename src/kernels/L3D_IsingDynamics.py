"""Kernel for Ising dynamics on 3D lattices.

Combines topology helpers from :mod:`kernels.L3D` with dynamics helpers
from :mod:`kernels.IsingDynamics`.  Engine-agnostic: the graph engine
(NX vs GT) is handled entirely by :func:`~kernels.L3D.prepare_lattice`.
"""

from __future__ import annotations

from typing import Any

from lrgsglib import *
from .L3D import get_l3d_out_suffix, prepare_lattice
from .IsingDynamics import (
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

__all__ = ["run_simulation"]

# Map observable prefix → extractor function
_EXTRACTORS = {
    "ene": extract_ene_data,
    "magn": extract_magn_data,
    "sout": extract_sout_data,
    "coeffs": extract_coeffs_data,
}


def run_simulation(args: Any) -> None:
    """Run L3D Ising dynamics with quenched + thermal averaging.

    Supports progressive checkpointing (``--save_frequency``) and
    resume from partial runs.  Each observable is saved to its own
    NPZ file with an ``_nt=<count>`` progress token.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments containing lattice, Ising, and averaging
        parameters.
    """
    n_quench = getattr(args, "number_of_averages", 1)
    n_thermal = getattr(args, "n_thermal", 1)
    save_freq = getattr(args, "save_frequency", 0)
    out_suffix = get_l3d_out_suffix(args)
    save_flags = resolve_save_flags(args)
    any_save = any(save_flags.values())

    # Single-quench array-task mode: if --quench-id >= 1, run only that
    # quench (for SLURM --array dispatch, one task per realisation).
    qid = getattr(args, "quench_id", -1)
    single_task_mode = qid >= 1
    if single_task_mode:
        quench_indices = [qid]
    else:
        quench_indices = list(range(1, n_quench + 1))

    # ── Pre-scan for completed quenches ──────────────────────────
    completed_q: set[int] = set()
    if any_save and n_quench > 1 and not single_task_mode:
        probe = prepare_lattice(args)
        # Stash topo_n_modes on lattice for stem builder
        rl = getattr(args, "runlang", "C1")
        if "topo" in rl.lower():
            probe._topo_n_modes = getattr(args, "topo_n_modes", 0)
        active_prefixes = [p for p, on in save_flags.items() if on]
        first_pfx = active_prefixes[0]
        stem0, save_dir0 = build_ising_stem(
            first_pfx, rl, probe, 1, is_l3d=True,
        )
        # stem without the _q=001 suffix
        stem_no_q = stem0.rsplit("_q=", 1)[0]
        completed_q = find_completed_quenches(
            save_dir0, stem_no_q, n_thermal,
            prefixes=active_prefixes,
        )
        if completed_q:
            print(
                f"Resuming: skipping {len(completed_q)} completed "
                f"quenches q={sorted(completed_q)}"
            )

    # ── Quench loop ──────────────────────────────────────────────
    for _q in quench_indices:
        if _q in completed_q:
            continue

        lattice = prepare_lattice(args)

        # Eigenvector-based init for ground_state_k conditions
        ic = getattr(args, "init_cond", "rand")
        if ic.startswith("ground_state") or ic.startswith("gs"):
            k = int(ic.rsplit("_", 1)[-1])
            if hasattr(lattice, "compute_k_eigvV") and hasattr(
                lattice, "load_eigV_on_graph"
            ):
                lattice.compute_k_eigvV(k + 1)
                lattice.load_eigV_on_graph(k, binarize=True)
            else:
                import warnings

                warnings.warn(
                    f"ground_state init requested but "
                    f"{type(lattice).__name__} does not support "
                    "eigenvector loading. Using random init.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                args.init_cond = "rand"

        # ── Set up per-observable accumulators ───────────────────
        rl = getattr(args, "runlang", "C1")
        if "topo" in rl.lower():
            lattice._topo_n_modes = getattr(args, "topo_n_modes", 0)

        obs_state: dict[str, dict[str, Any]] = {}
        for prefix, enabled in save_flags.items():
            if not enabled:
                continue
            stem, save_dir = build_ising_stem(
                prefix, rl, lattice, _q, is_l3d=True,
            )
            prior, n_done, old_path = load_ising_checkpoint(
                save_dir, stem,
            )
            acc = (
                rebuild_accumulator(prior)
                if prior is not None
                else new_accumulator(lattice, args, _q, is_l3d=True)
            )
            obs_state[prefix] = {
                "stem": stem,
                "dir": save_dir,
                "acc": acc,
                "n_done": n_done,
                "old_path": old_path,
            }

        # Resume point: minimum across all observables
        n_done = (
            min(s["n_done"] for s in obs_state.values())
            if obs_state
            else 0
        )

        # ── Thermal loop ────────────────────────────────────────
        for _t in range(n_done, n_thermal):
            isdy = run_ising_dynamics(
                args, lattice, out_suffix=out_suffix,
            )

            # Extract and accumulate each observable
            for prefix, state in obs_state.items():
                single = _EXTRACTORS[prefix](isdy)
                append_run(state["acc"], single)

            current_nt = _t + 1
            if should_checkpoint(current_nt, save_freq, n_thermal):
                for prefix, state in obs_state.items():
                    final = finalize_accumulator(
                        state["acc"], current_nt,
                    )
                    new_path = save_ising_checkpoint(
                        final,
                        state["dir"],
                        state["stem"],
                        current_nt,
                        old_path=state["old_path"],
                    )
                    state["old_path"] = new_path

            if args.remove_files:
                clean_up_files(isdy, lattice)
