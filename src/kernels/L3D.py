from typing import Any
from lrgsglib import Lattice3D, load_or_compute_Lattice3D

__all__ = [
        "initialize_l3d_dict_args",
        "prepare_lattice"]

def initialize_l3d_dict_args(args: Any) -> dict:
    """
    Generate a dictionary of arguments for initializing a Lattice3D object.

    This function processes the input arguments to create a dictionary
    containing the necessary parameters for initializing a Lattice3D instance.

    Parameters
    ----------
    args : Any
        Argument object containing lattice configuration, including dimensions,
        geometry, working directory, and flipping probability.

    Returns
    -------
    dict
        A dictionary with keys and values required for Lattice3D initialization.
    """
    argdict = dict(dim=args.L, geo=args.geometry, sgpathn=args.workdir, 
                    pflip=args.p)
    match args.cell_type:
        case 'rand':
            pass
        case _:
            argdict |= dict(init_nw_dict=True)
    return argdict

def prepare_lattice(args: Any, **kwargs) -> Any:
    """
    Initialize and configure a 3D lattice structure based on the provided 
    arguments. This function creates a Lattice3D object, applies modifications 
    to its edges, and optionally computes its Laplacian spectrum with 
    eigenvalues.

    Parameters
    ----------
    args : Any
        Argument object containing lattice configuration, including dimensions,
        geometry, and other properties.
    with_eigV : bool, optional
        If True, computes the Laplacian spectrum with eigenvalues for the
        lattice (default is True).

    Returns
    -------
    Any
        A configured Lattice3D instance with applied modifications.
    """
    lattice = load_or_compute_Lattice3D(
        **initialize_l3d_dict_args(args), 
        cell_type=args.cell_type,
        **kwargs
    )
    return lattice