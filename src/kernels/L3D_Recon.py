from lrgsglib import *
#
from .generic import *
from .L3D import *
from .IsingDynamics import *
from .Recon import *
#
def save_results(
    spinovp: List[np.ndarray],
    args: Any,
    save_dir: Path,
    loglogger: Any,
    syshapename: str = ''
) -> None:
    """
    Save spin overlap results to a compressed file.

    Parameters
    ----------
    spinovp : List[np.ndarray]
        List of spin overlap arrays from simulations.
    args : Any
        Argument object with attributes `mode`, `L`, `p`, and `T`.
    save_dir : Path
        Directory where the result file will be saved.
    loglogger : Any
        Logger instance for output.
    syshapename : str, optional (default='')
        Custom system shape name to include in the filename.

    Raises
    ------
    ValueError
        If `args.mode` is not supported.
    """
    # Construct the filename
    basename = f"smatch_{args.mode}"
    fname = '_'.join([basename, syshapename or args.L, f'p={args.p:.3g}', 
                      f'T={args.T:.3g}', args.out_suffix, f"avg={args.number_of_averages}"])
    save_path = save_dir / f"{fname}.npz"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    # Save results based on mode
    save_spin_overlap_results(spinovp, save_path, args, loglogger)


def run_single_recon_lattice(
    args: Any,
    loglogger: Any,
    spinovp: List[np.ndarray],
    max_factor: int = 2,
    start: int = 0,
    step: int = 1
) -> Path:
    """
    Run a single reconstruction simulation and append the results.

    Parameters
    ----------
    args : Any
        Argument object with simulation and lattice configuration.
    loglogger : Any
        Logger instance for status output.
    spinovp : List[np.ndarray]
        List to store spin overlap arrays from simulations.
    max_factor : int, optional
        Maximum factor for the spectral basis, by default 2.
    start : int, optional
        Starting index for the spectral basis, by default 0.
    step : int, optional
        Step size for the spectral basis, by default 1.

    Returns
    -------
    Path
        Path to the lattice output directory (lattice.path_lrgsg).
    """
    lattice = prepare_lattice(args)
    basis = lattice.get_sgspect_basis(max_factor, start, step)
    loglogger.debug("Laplacian spectrum and eigenbasis computed.")

    isdy = run_ising_dynamics(args, lattice)
    loglogger.debug("Ising dynamics run completed.")

    if args.remove_files:
        clean_up_files(isdy, lattice)

    spinovp_tmp = compute_spin_match_series(isdy.s, basis)
    loglogger.info(f"Reconstruction with first component: {spinovp_tmp[0]}")
    spinovp.append(spinovp_tmp)
    loglogger.debug("One reconstruction completed.")

    return lattice.path_lrgsg

def compute_reconstructions(
    args: Any,
    loglogger: Optional[Any] = None
) -> tuple[List[np.ndarray], Path]:
    """
    Perform multiple simulations and compute spin overlap reconstructions.

    Parameters
    ----------
    args : Any
        Argument object with simulation and lattice configuration.
    loglogger : Optional[Any], optional
        Logger instance for status output, by default None.

    Returns
    -------
    Tuple[List[np.ndarray], Path]
        A tuple containing:
        - List of spin overlap arrays from each simulation.
        - Path to the lattice output directory (lattice.path_lrgsg).
    """
    loglogger = setup_logger(loglogger)
    loglogger.info(f"Computing reconstructions for {args.number_of_averages} \
                   averages.")

    spinovp = []
    for _ in range(args.number_of_averages):
        lattice_path = run_single_recon_lattice(args, loglogger, spinovp)

    return spinovp, lattice_path



def compute_reconstructions_progressive(
    args: Any,
    save_dir: Path,
    loglogger: Optional[Any] = None
) -> Path:
    """
    Perform multiple simulations and compute spin overlap reconstructions,
    saving results progressively in batches.

    Parameters
    ----------
    args : Any
        Argument object with simulation and lattice configuration.
    save_dir : Path
        Directory where the result files will be saved.
    loglogger : Optional[Any], optional
        Logger instance for status output, by default None.

    Returns
    -------
    Path
        Path to the lattice output directory (lattice.path_lrgsg).
    """
    loglogger = setup_logger(loglogger)
    loglogger.info(f"Computing reconstructions for {args.number_of_averages} \
                   averages with progressive saving.")

    spinovp = []
    lattice = None  # will hold the last lattice object

    for i in range(args.number_of_averages):
        lattice_path = run_single_recon_lattice(args, loglogger, spinovp)
        # Save results incrementally to a single file, replacing the previous one
        if (i + 1) % args.freq == 0 or (i + 1) == args.number_of_averages:
            save_incremental(spinovp, (i + 1), args, save_dir, loglogger)

    return lattice_path 
