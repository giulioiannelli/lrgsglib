import pickle as pk
#
from numpy.random import randint
from numpy.typing import NDArray
from typing import Any, Union, Tuple
#
from ...config.funcs import peq_fstr
from .Lattice3D import Lattice3D, L3D_ONREP
#
__all__ = [
    "load_or_compute_Lattice3D"
]
#
def load_or_compute_Lattice3D(
        dim: Union[int, Tuple[int, int, int]], 
        geo: str = 'sc', 
        *,
        cell_type: str = 'rand',
        save: bool = True, 
        compute: str = 'energy_cupy',
        **kwargs: Any
) -> Lattice3D:
    """
    Load or compute a 3D lattice object, optionally saving it to a file.

    This function attempts to load a precomputed lattice object from a file.
    If the file does not exist, it computes the lattice, applies random edge
    flips, computes its energy and eigenvalues, and optionally saves the
    lattice object to a file.

    Parameters
    ----------
    side : int
        Linear size of the lattice (number of sites per edge).
    geo : str
        Geometry of the lattice (e.g., 'tri', 'hex', 'sqr').
    save : bool, optional
        If True, saves the computed lattice object to a file (default is True).
    **kwargs : Any
        Additional keyword arguments passed to the `Lattice2D` constructor.

    Returns
    -------
    Lattice2D
        A lattice object with computed energy and eigenvalues.

    Raises
    ------
    KeyError
        If the `pflip` argument is not provided in `kwargs`.

    Notes
    -----
    - The lattice object is saved to a file named based on its geometry and
      flipping probability (`pflip`).
    - The eigenvalues and energy are computed using the `cupy` backend.
    """
    tmp_l = Lattice3D(dim, geo, only_const_mode=True, **kwargs)
    if not compute:
        return tmp_l
    pflip = kwargs.get('pflip', .0)
    if cell_type != 'rand':
        kwargs.setdefault('init_nw_dict', True)
    new_seed = kwargs.get('seed', randint(0, 2**31 - 1)) if (pflip > 0. and pflip < 1.) else None
    seed_str = f'_seed{new_seed%2**16}' if new_seed is not None else ''
    match compute:
        case _ if compute.startswith('spectrum'):
            basename = 'L3D_spe'
        case _ if compute.startswith('energy'):
            basename = 'L3D_ene'
        case _ if compute.startswith('eigV'):
            howmany = int(compute.split('_')[1]) if '_' in compute else 1
            basename = 'L3D_eigV_' + str(howmany)
        case _:
            raise ValueError(f"Unknown compute option '{compute}'.")
    fname = '_'.join([basename, peq_fstr(pflip)]) + seed_str
    #
    pname = tmp_l.path_graph / (fname + '.pkl')
    #
    if pname.exists():
        lattice = pk.load(open(pname, 'rb'))
        lattice.__init_loaded_graph__(path_data=kwargs.get('path_data', None))
    else:
        lattice = Lattice3D(dim, geo=geo, **kwargs)
        match cell_type:
            case 'rand':
                lattice.flip_random_fract_edges()
            case _:
                try:
                    pattern = lattice.nwDict[cell_type][L3D_ONREP]
                except KeyError:
                    raise ValueError(f"Unknown cell_type '{cell_type}'.") from None
                lattice.flip_sel_edges(pattern)
        try:
            import cupy as cp
            # Test if a GPU is available
            cp.cuda.runtime.getDeviceCount()
        except Exception as e:
            compute = compute.replace('cupy', 'numpy')
        match compute:
            case _ if compute.startswith('spectrum'):
                routine = compute.split('_')[-1] if '_' in compute else 'cupy'
                lattice.compute_laplacian_spectrum_weigV(with_routine=routine)
            case _ if compute.startswith('energy'):
                routine = compute.split('_')[-1] if '_' in compute else 'cupy'
                lattice.compute_rbim_energy_eigV_all(with_routine=routine)
            case _ if compute.startswith('eigV'):
                routine = compute.split('_')[-1] if '_' in compute else 'cupy'
                lattice.compute_k_eigvV(k=howmany, with_routine=routine)
            case _:
                raise ValueError(f"Unknown compute option '{compute}'.")
        if save:
            with open(pname, 'wb') as f:
                pk.dump(lattice, f)
    return lattice