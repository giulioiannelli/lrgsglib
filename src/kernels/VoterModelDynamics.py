"""VoterModel dynamics runner functions.

Light runner around :class:`lrgsglib.statsys.VoterModel`: the model owns its own
observable persistence (``savedisk=True`` writes ``m``/``sout``/``cldist`` under
the lattice ``dynpath``), so -- unlike IsingDynamics -- there is no NPZ
accumulator / checkpoint suite here. Each call runs ONE realisation and lets the
model write its files; the quench index is folded into ``out_suffix`` so
independent realisations land in distinct files.

Command example (2D lattice, linear voter, Python backend, 200 realisations)::

    python lrgsglib/src/L2D_VoterModel.py 64 -p 0.1 -na 200 -wd vmtest \\
        --rule linear --upd_mode asynchronous -rl py --simref 1000 --save-magn
"""

from typing import Any

from lrgsglib import join_non_empty
from lrgsglib.statsys import VoterModel

__all__ = [
    "get_out_suffix",
    "initialize_voter_dict_args",
    "run_voter_model",
]


def get_out_suffix(args: Any, quench_id: int | None = None) -> str:
    """Build the per-run output suffix (folds in the quench index).

    The voter model keys its output filenames on ``out_suffix``; folding the
    quench index in (``q<id>``) keeps independent realisations from overwriting
    one another and makes the array (``-qid``) and in-process loop modes produce
    identically named files for the same realisation index.
    """
    base = join_non_empty("_", args.init_cond, args.cell_type, args.out_suffix)
    if quench_id is not None and quench_id >= 0:
        return join_non_empty("_", base, f"q{quench_id}")
    return base


def initialize_voter_dict_args(args: Any, out_suffix: str) -> dict[str, Any]:
    """Map parsed CLI args onto the VoterModel constructor kwargs.

    ``--save-all`` is a convenience switch that turns on every observable gate
    (magnetization, trajectory, cluster-size distribution).
    """
    save_all = bool(getattr(args, "save_all", False))
    return dict(
        # Axis A -- local update rule
        rule=args.rule,
        q=args.q,
        eps=args.eps,
        alpha=args.alpha,
        # Axis B -- update schedule
        upd_mode=args.upd_mode,
        # backend / init / time
        runlang=args.runlang,
        ic=args.init_cond,
        steps=args.steps,
        simref=args.simref,
        state_type=args.state_type,
        # absorbing-state detection
        absorbing_check=args.absorbing_check,
        absorbing_every=args.absorbing_every,
        # cluster-size distribution
        track_clusters=save_all or bool(args.save_cldist),
        cluster_mode=args.cluster_mode,
        # observable gates + trajectory subsampling
        savemagn=save_all or bool(args.save_magn),
        savedyn=save_all or bool(args.save_sout),
        sout_every=args.sout_every,
        sout_nlog=args.sout_nlog,
        sout_force_full=args.sout_force_full,
        # sampling / output
        freq=args.freq,
        nSampleLog=args.num_log_samples,
        seed=args.seed,
        rndStr=args.randstr,
        out_suffix=out_suffix,
    )


def run_voter_model(
    args: Any, signed_graph: Any, quench_id: int | None = None,
) -> VoterModel:
    """Instantiate and run one VoterModel realisation on ``signed_graph``.

    The model persists its enabled observables to disk itself (``savedisk``);
    the C backend's transient export files are removed by ``run(clean_export=...)``
    when ``--remove_files`` is set.
    """
    out_suffix = get_out_suffix(args, quench_id)
    vm = VoterModel(signed_graph, **initialize_voter_dict_args(args, out_suffix))
    vm.run(
        tqdm_on=False,
        verbose=getattr(args, "verbose", False),
        clean_export=bool(getattr(args, "remove_files", True)),
    )
    return vm
