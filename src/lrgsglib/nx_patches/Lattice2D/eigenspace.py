import pickle as pk
#
from numpy.random import randint
from numpy.typing import NDArray
from typing import Any
#
from ...config.funcs import peq_fstr
from .Lattice2D import Lattice2D, L2D_ONREP
#
__all__ = [
    "create_lattice_with_eigenspace",
    "load_or_compute_Lattice2D"
]
#
def create_lattice_with_eigenspace(
    side: int,
    disorder_struct: str = 'random',
    *,
    with_routine: str = 'scipy',
    k: int = 1,
    **kwargs: Any
) -> Lattice2D:
    """
    Construct a 2D lattice, apply a specified disorder pattern, and compute its first k eigenpairs.

    Parameters
    ----------
    side : int
        Linear size of the lattice (number of sites per edge).
    disorder_struct : {'random'} ∪ str, default 'random'
        Name of the disorder structure to apply. If 'random', a random fraction
        of edges will be flipped. Otherwise, must be a key in
        `l.nwDict[disorder_struct]` to select a predefined pattern.
    with_routine : str, keyword-only, default 'scipy'
        Backend routine to compute eigenvalues/vectors (e.g., 'scipy', 'numpy').
    k : int, keyword-only, default 1
        Number of lowest-magnitude eigenvalues (and corresponding eigenvectors)
        to compute.
    **kwargs : Any, keyword-only
        Additional keyword arguments passed to `Lattice2D` constructor. If
        `disorder_struct != 'random'`, the flag `init_nw_dict=True` will be
        injected to enable loading from `nwDict`.

    Returns
    -------
    Lattice2D
        An initialized lattice object with `.eigv` and `.eigV` attributes
        holding the computed eigenvalues and eigenvectors.

    Raises
    ------
    ValueError
        If `disorder_struct` is not 'random' and not found in the lattice’s
        `nwDict`.
    """
    # Ensure that predefined structures load their network definitions
    if disorder_struct != 'random':
        kwargs.setdefault('init_nw_dict', True)

    # Initialize lattice
    lattice = Lattice2D(side, **kwargs)

    # Apply disorder
    if disorder_struct == 'random':
        lattice.flip_random_fract_edges()
    else:
        try:
            pattern = lattice.nwDict[disorder_struct][L2D_ONREP]
        except KeyError:
            raise ValueError(f"Unknown disorder_struct '{disorder_struct}'.") from None
        lattice.flip_sel_edges(pattern)

    # Compute the lowest-k eigenpairs
    lattice.compute_k_eigvV(with_routine=with_routine, k=k)
    return lattice

def load_or_compute_Lattice2D(
        side1: int, 
        geo: str,
        *,
        cell_type: str = 'rand',
        save: bool = True, 
        compute: str = 'energy_cupy',
        **kwargs: Any
) -> Lattice2D:
    """
    Load or compute a 2D lattice object, optionally saving it to a file.

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
    compute : str, optional
        Specifies the computation to perform on the lattice. Options include:
        - 'spectrum_cupy': Compute the spectrum using cupy.
        - 'energy_cupy': Compute the energy using cupy.
        - 'spectrum_scipy': Compute the spectrum using scipy.
        - 'energy_scipy': Compute the energy using scipy.
        Default is 'energy_cupy'.
    cell_type : str, optional
        Type of cell structure to use for the lattice. If 'rand', a random
        fraction of edges will be flipped. Otherwise, it must be a key in
        `lattice.nwDict[cell_type]` to select a predefined pattern.
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
    tmp_l = Lattice2D(side1, geo, only_const_mode=True, **kwargs)
    if not compute:
        return tmp_l
    pflip = kwargs.get('pflip', .0)
    new_seed = kwargs.get('seed', randint(0, 2**31 - 1)) if (pflip > 0. and pflip < 1.) else None
    seed_str = f'_seed{new_seed%2**16}' if new_seed is not None else ''
    match compute:
        case _ if compute.startswith('spectrum'):
            basename = 'L2D_spe'
        case _ if compute.startswith('energy'):
            basename = 'L2D_ene'
    fname = '_'.join([basename, peq_fstr(pflip)]) + seed_str
    #
    pname = tmp_l.path_graph / (fname + '.pkl')
    #
    if pname.exists():
        lattice = pk.load(open(pname, 'rb'))
        lattice.__init_loaded_graph__(path_data=kwargs.get('path_data', None))
    else:
        lattice = Lattice2D(side1, geo=geo, **kwargs)
        match cell_type:
            case 'rand':
                lattice.flip_random_fract_edges()
            case _:
                try:
                    pattern = lattice.nwDict[cell_type][L2D_ONREP]
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
                routine = compute.split('_')[1] if '_' in compute else 'numpy'
                lattice.compute_laplacian_spectrum_weigV(with_routine=routine)
            case _ if compute.startswith('energy'):
                routine = compute.split('_')[1] if '_' in compute else 'numpy'
                lattice.compute_rbim_energy_eigV_all(with_routine=routine)
            case _:
                raise ValueError(f"Unknown compute option '{compute}'.")
        if save:
            with open(pname, 'wb') as f:
                pk.dump(lattice, f)
    return lattice