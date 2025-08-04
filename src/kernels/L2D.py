from typing import Any
from lrgsglib import Lattice2D, load_or_compute_Lattice2D

__all__ = [
    "initialize_l2d_dict_args",
    "prepare_lattice"
]

def initialize_l2d_dict_args(args: Any) -> dict:
    """
    Generate a dictionary of arguments for initializing a Lattice2D object.

    Parameters:
    -----------
    args : Any
        Argument object containing lattice configuration.

    Returns:
    --------
    dict
        Dictionary of arguments for Lattice2D initialization.

    Notes:
    ------
    The function determines whether to include the `init_nw_dict` key based on the `cell_type` attribute of `args`. If `cell_type` is 'rand', the dictionary excludes `init_nw_dict`. Otherwise, it includes `init_nw_dict` set to `True`.
    """
    match args.cell_type:
        case 'rand':
            return dict(side1=args.L, geo=args.geometry, sgpathn=args.workdir, 
                    pflip=args.p)
        case _:
            return dict(side1=args.L, geo=args.geometry, sgpathn=args.workdir, 
                    pflip=args.p, init_nw_dict=True)
        
def prepare_lattice(args: Any, **kwargs) -> Any:
    """
    Initialize and modify a Lattice2D object based on args.

    Parameters
    ----------
    args : Any
        Argument object containing lattice configuration.

    Returns
    -------
    Any
        Configured Lattice2D instance.
    """
    lattice = load_or_compute_Lattice2D(
        **initialize_l2d_dict_args(args), 
        cell_type=args.cell_type,
        **kwargs
    )
    return lattice

