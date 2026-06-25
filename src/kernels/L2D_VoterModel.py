"""2D-lattice driver for the VoterModel program.

Builds a fresh disorder realisation per quench (a direct ``Lattice2D`` +
``flip_random_fract_edges`` -- NOT the cached ``load_or_compute`` path, so each
realisation has independent disorder) and runs the voter model on it. Two
dispatch modes:

* ``--quench-id N`` (>= 0): run exactly realisation ``N`` (one SLURM array task).
* otherwise (-1, default): loop over ``--number_of_averages`` realisations.

The model persists its own observables; this driver only owns lattice
construction and the realisation loop.
"""

from lrgsglib import Lattice2D

from .L2D import initialize_l2d_dict_args, prepare_lattice as _prepare_lattice_gt
from .VoterModelDynamics import run_voter_model
from parsers.shared import get_graph_engine

__all__ = ["run_simulation"]


def _prepare_lattice(args):
    """Return a fresh 2D signed lattice (independent disorder each call)."""
    engine = get_graph_engine(args)
    if engine == "gt":
        return _prepare_lattice_gt(args)

    # Default NX path: direct construction + a new random flip realisation
    # (the cached load_or_compute path would reuse one disorder across quenches).
    lattice = Lattice2D(**initialize_l2d_dict_args(args))
    if lattice.init_nw_dict:
        lattice.flip_sel_edges(lattice.nwDict[args.cell_type]["G"])
    else:
        lattice.flip_random_fract_edges()
    return lattice


def _run_one(args, quench_id):
    lattice = _prepare_lattice(args)
    return run_voter_model(args, lattice, quench_id=quench_id)


def run_simulation(args):
    qid = getattr(args, "quench_id", -1)
    if qid is not None and qid >= 0:
        # Single realisation -- one SLURM array task.
        _run_one(args, qid)
        return
    # In-process loop over independent disorder realisations.
    for i in range(1, args.number_of_averages + 1):
        _run_one(args, i)
