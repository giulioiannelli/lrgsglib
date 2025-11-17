from lrgsglib import Lattice2D

from .ContactProcessDynamics import run_contact_process
from .L2D import initialize_l2d_dict_args


def _prepare_lattice(args):
    lattice = Lattice2D(**initialize_l2d_dict_args(args))
    if lattice.init_nw_dict:
        lattice.flip_sel_edges(lattice.nwDict[args.cell_type]["G"])
    else:
        lattice.flip_random_fract_edges()
    return lattice


def run_simulation(args):
    for _ in range(args.number_of_averages):
        lattice = _prepare_lattice(args)
        run_contact_process(args, lattice)
