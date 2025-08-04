import numpy as np
#
from typing import Any, List, Callable, Optional, Union, Tuple
from pathlib import Path
#
from numpy.typing import NDArray
#
from lrgsglib.config.funcs import avgeq_fstr, build_fname_or_pattern_direct
from lrgsglib.utils.basic import update_mean_m2, join_non_empty
from lrgsglib.utils.lrg import compute_spin_match_series
#
from .generic import setup_logger
from .IsingDynamics import run_ising_dynamics, clean_up_files
#
__all__ = [
    "build_fname_or_pattern",
    "check_existing_file_and_needed_averages",
    "run_reconstruction",
    "save_spin_overlap_results",
    "save_incremental",
    "compute_recon_prog_lattice_incr",
]

recon_fname = 'smatch'
#
def build_fname_or_pattern(
    args: Any,
    mode: str = 'fname',
    n_avg: int = 0,
    ext: str = '.npz',
    basefn: str = recon_fname
) -> str:
    """
    Build a filename or pattern for saving results based on the provided 
    arguments.

    Parameters
    ----------
    args : Any
        Argument object with attributes `mode`, `L`, `p`, `T`, and `out_suffix`.
    mode : str, optional
        Mode for the filename, either 'fname' or 'pattern', by default 'fname'.
    n_avg : int, optional
        Number of averages to include in the filename, by default 0.
    with_ext : bool, optional
        Whether to include the file extension in the filename, by default True.
    Returns
    -------
    str
        The constructed filename or pattern.
    """
    basename = '_'.join([basefn, args.mode])
    return build_fname_or_pattern_direct(basename, args.p, args.T, args.out_suffix, mode, n_avg, ext)
#


#
def check_existing_file_and_needed_averages(
    save_dir: Path,
    args: Any,
    loglogger: Any
) -> int:
    """
    Check if the required file exists and determine the number of additional
    averages needed to reach the requested total.

    Parameters
    ----------
    save_dir : Path
        Directory where the result file is expected to be saved.
    args : Any
        Argument object with attributes `mode`, `L`, `p`, `T`, and 
        `number_of_averages`.
    loglogger : Any
        Logger instance for status output.

    Returns
    -------
    int
        The number of additional averages needed to reach the requested total.
        Returns 0 if the file already satisfies the requested # of averages.
    """
    pattern = build_fname_or_pattern(args, mode='pattern')
    existing_files = list(save_dir.glob(pattern))
    #
    if not existing_files:
        loglogger.info(f"No existing file found. Computing \
                       {args.number_of_averages} averages.")
        return args.number_of_averages
    #
    # Extract the current number of averages from the filename
    existing_file = existing_files[0]
    loglogger.info(f"Found existing file: {existing_file}")
    try:
        current_avg = int(existing_file.stem.split("_avg=")[-1])
    except ValueError:
        loglogger.warning("Could not parse the number of averages from the \
                          filename.")
        return args.number_of_averages
    #
    # Determine the number of additional averages needed
    if current_avg >= args.number_of_averages:
        loglogger.info("The existing file already satisfies the requested \
                       averages.")
        return 0
    needed_averages = args.number_of_averages - current_avg
    loglogger.info(f"Existing file has {current_avg} averages. "
                   f"{needed_averages} more averages are needed.")
    return needed_averages
#
def run_reconstruction(
    args: Any,
    loglogger: Any,
    spinovp: List[NDArray],
    prepare_lattice_func: Callable[[Any], Any],
) -> Union[List[NDArray], Path]:
    """
    Run a reconstruction simulation with customizable body functions.

    Parameters
    ----------
    args : Any
        Argument object with simulation and lattice configuration.
    loglogger : Any
        Logger instance for status output.
    spinovp : List[np.ndarray]
        List to store spin overlap arrays from simulations.
    prepare_lattice_func : Callable[[Any], Any]
        Function to prepare the lattice.
    max_factor : int, optional
        Maximum factor for the basis computation, by default 2.

    Returns
    -------
    Path
        Path to the lattice output directory (lattice.path_lrgsg).
    """
    pflip = args.p
    assert pflip >= 0. and pflip <= 1., "pflip must be in [0, 1]"
    if pflip > 0. and pflip < 1.:
        lattice = prepare_lattice_func(args, save=False)
    else:
        lattice = prepare_lattice_func(args)
    
    basis = lattice.get_sgspect_basis(args.max_factor, args.start_idx_basis, args.basis_step)
    loglogger.debug("Basis computed.")

    isdy = run_ising_dynamics(args, lattice)
    loglogger.debug("Dynamics run completed.")
    if args.remove_files:
        clean_up_files(isdy, lattice)

    spinovp_tmp = compute_spin_match_series(isdy.s, basis)
    loglogger.info(f"Recon. w. 1st component: {spinovp_tmp[0]}")
    match args.mode:
        case 'statistic':
            if not spinovp:
                basis_len = len(spinovp_tmp)
                spinovp = [np.zeros(basis_len), np.zeros(basis_len), 0]
            spinovp = update_mean_m2(*spinovp, spinovp_tmp)
        case 'full': 
            spinovp.append(spinovp_tmp)
    loglogger.debug("Reconstruction completed.")

    return spinovp, lattice.path_lrgsg
#
def save_spin_overlap_results(
    spinovp: List[NDArray],
    save_path: Path,
    args: Any,
    loglogger: Any
) -> None:
    """
    Save spin overlap results to a compressed file.

    Parameters
    ----------
    spinovp : List[np.ndarray]
        List of spin overlap arrays from simulations.
    save_path : Path
        Path where the result file will be saved.
    args : Any
        Argument object with attributes `mode`, `L`, `p`, and `T`.
    loglogger : Any
        Logger instance for output.

    Raises
    ------
    ValueError
        If `args.mode` is not supported.
    """

    # Save results based on mode
    if args.mode == 'statistic':
        np.savez_compressed(save_path, mean=spinovp[0], std=np.sqrt(spinovp[1]/spinovp[2]))
    elif args.mode == 'full':
        data = np.stack(spinovp)
        np.savez_compressed(save_path, sovstack=data)
    else:
        raise ValueError(f"Unsupported mode '{args.mode}'.")

    loglogger.info(f"Saved result to: {save_path}")

def save_incremental(
    spinovp: List[NDArray],
    current_avg: int,
    args: Any,
    save_dir: Path,
    loglogger: Any
) -> None:
    """
    Save spin overlap results incrementally to a single file, replacing the 
    previous file.

    Parameters
    ----------
    spinovp : List[np.ndarray]
        List of spin overlap arrays from simulations.
    current_avg : int
        The current number of averages computed.
    args : Any
        Argument object with attributes `mode`, `L`, `p`, and `T`.
    save_dir : Path
        Directory where the result file will be saved.
    loglogger : Any
        Logger instance for output.
    """
    n_avg = current_avg - args.save_frequency
    fname = build_fname_or_pattern(args, n_avg=n_avg)
    save_path = save_dir / fname
    loglogger.debug(f"Preparing to save results to: {save_path}")
    # Ensure the save directory exists
    # Remove the previously saved file if it exists
    if save_path.exists():
        loglogger.debug(f"Removing existing file: {save_path}")
        save_path.unlink()
    n_avg_new = current_avg
    fname = build_fname_or_pattern(args, n_avg=n_avg_new)
    new_save_path = save_dir / fname
    # Save the new results
    save_spin_overlap_results(spinovp, new_save_path, args, loglogger)

def compute_recon_prog_lattice_incr(
    args: Any,
    loglogger: Optional[Any],
    prepare_lattice_func: Callable[[Any], Any],
) -> Path:
    """
    Perform multiple simulations and compute spin overlap reconstructions,
    saving results progressively in batches. Only runs the needed number of
    simulations based on existing results.

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
    loglogger.info(f"Checking existing results for {args.number_of_averages} \
                   averages.")
    # Check existing file and determine needed averages
    test_lattice = prepare_lattice_func(args, compute = None)
    save_dir = test_lattice.path_lrgsg  # Prepare lattice without eigV
    needed_averages = check_existing_file_and_needed_averages(save_dir, args, loglogger)

    spinovp = []
    if needed_averages == 0:
        loglogger.info("No additional simulations needed. Exiting.")
        return save_dir
    elif needed_averages < args.number_of_averages:
        loglogger.info(f"Only {needed_averages} additional simulations needed.")
        given_averages = args.number_of_averages - needed_averages
        fname_read = build_fname_or_pattern(args, n_avg=given_averages)
        loglogger.debug(f"Loading existing results from: {fname_read}")
        path_read = save_dir / fname_read
        data = np.load(path_read)
        # Load existing results if available
        match args.mode:
            case 'statistic':
                spinovp = [data['mean'], data['std']**2*given_averages, given_averages]
            case 'full':
                spinovp = list(data['sovstack'])      # Load full stack
        loglogger.info(f"Loaded {given_averages} existing averages from: {path_read}")
        given_averages -= args.save_frequency  # Adjust given averages to account for the next batch
    elif needed_averages == args.number_of_averages:
        given_averages = 0
        loglogger.debug("No existing results found. Starting from scratch.")
    else:
        raise ValueError(f"Unexpected number of needed averages: {needed_averages}")
    # Perform additional simulations
    loglogger.info(f"Running {needed_averages} additional simulations.")
    lattice_path = None  # will hold the last lattice path
    for i in range(needed_averages):
        spinovp, lattice_path = run_reconstruction(args, loglogger, spinovp, prepare_lattice_func)
        # Save results incrementally to a single file, replacing the previous one
        given_averages += 1
        if given_averages % args.save_frequency == 0:
            save_incremental(spinovp, given_averages, args, save_dir, loglogger)
    save_incremental(spinovp, args.number_of_averages, args, save_dir, loglogger)  # Remove the file after loading

    return lattice_path