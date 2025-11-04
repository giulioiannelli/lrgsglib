from typing import Any
from numpy.typing import NDArray
from lrgsglib import Lattice2D, load_or_compute_Lattice2D, flip_to_positive_majority

__all__ = [
    "initialize_l2d_dict_args",
    "prepare_lattice",
    "eigV_for_lattice2D",
    "adjust_eigV_for_lattice2D",
    "eigV_for_lattice2D_ptch",
    "eigv_for_lattice2D",
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




def eigV_for_lattice2D(side, mode='scipy', howmany=1, **kwargs) -> NDArray:
    l = Lattice2D(side, **kwargs)
    l.flip_random_fract_edges()
    l.compute_k_eigvV(backend=mode, k=howmany)
    return l.eigV

def adjust_eigV_for_lattice2D(leigV: NDArray) -> NDArray:
    for i,leV in enumerate(leigV):
        leigV[i] = flip_to_positive_majority(leV)
    return leigV

def eigV_for_lattice2D_ptch(**kwargs) -> NDArray:
    return adjust_eigV_for_lattice2D(eigV_for_lattice2D(**kwargs))

def eigv_for_lattice2D(side, mode: str = "full", **kwargs) -> NDArray:
    l = Lattice2D(side, **kwargs)
    l.flip_random_fract_edges()
    match mode:
        case "full":
            l.compute_laplacian_spectrum_weigV()
        case _ if mode.startswith("some"):
            k = int(mode.split("_")[-1])
            l.compute_k_eigvV(k=k)
    return l.eigv

