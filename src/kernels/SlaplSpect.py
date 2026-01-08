"""
Generic spectral analysis framework for signed Laplacian computations.

This module provides reusable utilities for computing and analyzing eigenvalue
and eigenvector distributions across different graph types (Lattice2D, MCG, etc.).

The framework uses callback functions to handle graph-specific operations while
providing generic logic for:
- Incremental computation with checkpointing
- Distribution binning and statistics
- File I/O and data persistence
- Resume from partial results
"""

import numpy as np
import pickle as pk
from pathlib import Path
from collections import Counter
from typing import Any, Callable, Optional, List, Tuple, Union

from lrgsglib.utils.basic.probability import (
    create_symmetric_log_bins,
    linear_binning_hist
)
from lrgsglib.config.funcs import bin_eigenvalues
from .generic import setup_logger

__all__ = [
    "process_eigen_distribution",
    "save_data",
    "load_existing_data",
    "check_existing_and_needed_averages",
    "create_initial_bin_counter",
    "build_eigval_fname_base",
    "build_eigvec_fname_base",
    "find_existing_data_by_na",
    "save_with_na",
    "cleanup_intermediate_files",
]


def save_data(args: Any, data: Any, fname: Path) -> None:
    """
    Save computed data to a pickle file.

    Parameters
    ----------
    args : Any
        Argument object with `verbose` attribute for logging.
    data : Any
        Data to save (typically Counter or list of arrays).
    fname : Path
        Full path to the output file.

    Returns
    -------
    None
        Data is saved to disk.

    Examples
    --------
    >>> from collections import Counter
    >>> save_data(args, Counter([1, 2, 3]), Path("output.pkl"))
    """
    if args.verbose:
        print(f"Saving data to {fname}...")
    with open(fname, "wb") as f:
        pk.dump(data, f)
    if args.verbose:
        print(f"Data saved to {fname}")


def load_existing_data(
    path_fname: Path,
    mode: str,
    verbose: bool = False
) -> Tuple[Any, np.ndarray]:
    """
    Load existing checkpoint data from file.

    Parameters
    ----------
    path_fname : Path
        Path to the checkpoint file.
    mode : str
        Processing mode: "eigvec_dist" or "eigval_dist".
    verbose : bool, optional
        Enable logging output.

    Returns
    -------
    bin_counter : Any
        Loaded bin counter (Counter or list of Counter).
    initial_data : np.ndarray
        Reconstructed initial data for binning as numpy array.

    Examples
    --------
    >>> counter, data = load_existing_data(
    ...     Path("checkpoint.pkl"), "eigval_dist"
    ... )
    """
    with open(path_fname, "rb") as f:
        bin_counter = pk.load(f)

    if mode == "eigvec_dist":
        # Reconstruct flattened data from list of counters
        initial_data = np.array([val for counter in bin_counter for val in counter.elements()])
    else:
        # Reconstruct from single counter
        initial_data = np.array(list(bin_counter.elements()))

    if verbose:
        print(f"Loaded existing data from {path_fname}")

    return bin_counter, initial_data


def check_existing_and_needed_averages(
    working_path: Path,
    fname_base: str,
    verbose: bool = False
) -> Tuple[int, Path]:
    """
    Check for existing checkpoint files and determine progress.

    Parameters
    ----------
    working_path : Path
        Directory to search for checkpoints.
    fname_base : str
        Base filename pattern (e.g., "mc_dist_eigval_p=0.1").
    verbose : bool, optional
        Enable logging output.

    Returns
    -------
    navg_done : int
        Number of averages already completed (0 if none found).
    path_fname : Path
        Path to existing checkpoint or new file location.

    Examples
    --------
    >>> navg, path = check_existing_and_needed_averages(
    ...     Path("data/"), "mc_dist_eigval_p=0.1"
    ... )
    """
    search_pattern = f"{fname_base}_na=*.pkl"
    existing_files = sorted(working_path.glob(search_pattern))

    navg_done = 0
    if existing_files:
        navg_done = max(
            int(f.stem.split('_')[-1].split('=', 1)[1])
            for f in existing_files
        )

    path_fname = working_path / f"{fname_base}_na={navg_done}.pkl"

    if verbose and navg_done > 0:
        print(f"Found existing checkpoint with {navg_done} averages")

    return navg_done, path_fname


def create_initial_bin_counter(mode: str, howmany: int = 0) -> Union[Counter, List[Counter]]:
    """
    Create initial bin counter structure based on mode.

    Parameters
    ----------
    mode : str
        Processing mode: "eigvec_dist" or "eigval_dist".
    howmany : int, optional
        Number of eigenvectors (for eigvec_dist mode).

    Returns
    -------
    bin_counter : Counter or list of Counter
        Empty counter(s) for accumulating binned values.

    Examples
    --------
    >>> counter = create_initial_bin_counter("eigval_dist")
    >>> counters = create_initial_bin_counter("eigvec_dist", howmany=5)
    """
    if mode == "eigvec_dist":
        return [Counter() for _ in range(howmany)]
    else:
        return Counter()


def build_eigval_fname_base(
    graph_prefix: str,
    pflip: float,
    cell_type: str = "",
    **kwargs
) -> str:
    """
    Build base filename for eigenvalue distribution results.

    Parameters
    ----------
    graph_prefix : str
        Graph type prefix (e.g., "", "mc_", "l3d_").
    pflip : float
        Edge flip probability.
    cell_type : str, optional
        Cell type identifier (for lattices).
    **kwargs
        Additional parameters (ignored for compatibility).

    Returns
    -------
    fname_base : str
        Base filename without extension.

    Examples
    --------
    >>> build_eigval_fname_base("mc_", 0.1)
    'mc_dist_eigval_p=0.1'
    >>> build_eigval_fname_base("", 0.05, cell_type="square")
    'dist_eigval_p=0.05_square'
    """
    if graph_prefix:
        # MCG style: "mc_dist_eigval_p=0.1"
        return f"{graph_prefix}dist_eigval_p={pflip:.3g}"
    else:
        # Lattice style: "dist_eigval_p=0.05_square"
        if cell_type:
            return f"dist_eigval_p={pflip:.3g}_{cell_type}"
        else:
            return f"dist_eigval_p={pflip:.3g}"


def build_eigvec_fname_base(
    graph_prefix: str,
    howmany: int,
    pflip: float,
    eigen_mode: str,
    **kwargs
) -> str:
    """
    Build base filename for eigenvector distribution results.

    Parameters
    ----------
    graph_prefix : str
        Graph type prefix (e.g., "", "mc_", "l3d_").
    howmany : int
        Number of eigenvectors analyzed.
    pflip : float
        Edge flip probability.
    eigen_mode : str
        Eigenvector selection mode ("smallest", "largest", "middle").
    **kwargs
        Additional parameters (ignored for compatibility).

    Returns
    -------
    fname_base : str
        Base filename without extension.

    Examples
    --------
    >>> build_eigvec_fname_base("mc_", 5, 0.1, "smallest")
    'mc_dist5_p=0.1_smallest'
    >>> build_eigvec_fname_base("", 10, 0.05, "largest")
    'dist10_p=0.05_largest'
    """
    if graph_prefix:
        # MCG style: "mc_dist5_p=0.1_smallest"
        return f"{graph_prefix}dist{howmany}_p={pflip:.3g}_{eigen_mode}"
    else:
        # Lattice style: "dist10_p=0.05_largest"
        return f"dist{howmany}_p={pflip:.3g}_{eigen_mode}"


def find_existing_data_by_na(
    working_path: Path,
    fname_base: str,
    requested_na: int,
    verbose: bool = False
) -> Tuple[Optional[Path], int, bool]:
    """
    Find existing data files using glob pattern and na values.

    Scans for files matching `{fname_base}_na=*.pkl` and determines:
    - If enough data already exists (na >= requested)
    - Which file to resume from (highest na < requested)

    Parameters
    ----------
    working_path : Path
        Directory to search for existing files.
    fname_base : str
        Base filename pattern (e.g., "mc_eigvals_p=0.0").
    requested_na : int
        Number of averages requested.
    verbose : bool, optional
        Enable logging output.

    Returns
    -------
    existing_file : Path or None
        Path to file to resume from (highest na < requested), or None if starting fresh.
    existing_na : int
        Number of averages in existing file (0 if none).
    should_skip : bool
        True if file with na >= requested already exists (computation not needed).

    Examples
    --------
    >>> file, na, skip = find_existing_data_by_na(
    ...     Path("data/"), "mc_eigvals_p=0.0", 20, verbose=True
    ... )
    Found existing data with 10 realizations
    Resuming from realization 10 to reach 20 total
    """
    glob_pattern = f"{fname_base}_na=*.pkl"
    matching_files = list(working_path.glob(glob_pattern))

    if not matching_files:
        if verbose:
            print(f"No existing data found. Computing {requested_na} realizations from scratch")
        return None, 0, False

    files_with_enough = []
    files_to_resume = []

    for fpath in matching_files:
        try:
            fname = fpath.stem
            na_part = fname.split('na=')[1]
            na_val = int(na_part)
            if na_val >= requested_na:
                files_with_enough.append((fpath, na_val))
            else:
                files_to_resume.append((fpath, na_val))
        except (IndexError, ValueError):
            continue

    # Check if we already have enough data
    if files_with_enough:
        files_with_enough.sort(key=lambda x: x[1])
        existing_file, existing_na = files_with_enough[0]
        if verbose:
            print(f"Found existing data with {existing_na} realizations at {existing_file}")
            print(f"Already have >= {requested_na} requested realizations. Skipping computation.")
        return existing_file, existing_na, True

    # Find highest na < requested to resume from
    if files_to_resume:
        files_to_resume.sort(key=lambda x: x[1], reverse=True)
        existing_file, existing_na = files_to_resume[0]
        if verbose:
            print(f"Found existing data with {existing_na} realizations at {existing_file}")
            print(f"Resuming from realization {existing_na} to reach {requested_na} total")
        return existing_file, existing_na, False

    if verbose:
        print(f"No usable existing data found. Computing {requested_na} realizations from scratch")
    return None, 0, False


def save_with_na(
    data: Any,
    working_path: Path,
    fname_base: str,
    current_na: int,
    metadata: Any = None,
    verbose: bool = False
) -> Path:
    """
    Save data with na suffix in filename.

    Parameters
    ----------
    data : Any
        Data to save (typically list of arrays or Counter objects).
    working_path : Path
        Directory to save to.
    fname_base : str
        Base filename (e.g., "mc_eigvals_p=0.0").
    current_na : int
        Number of averages in this save.
    metadata : Any, optional
        Optional metadata to save in separate file.
    verbose : bool, optional
        Enable logging output.

    Returns
    -------
    fname : Path
        Path to saved data file.

    Examples
    --------
    >>> save_with_na(
    ...     eigvals_list, Path("data/"), "mc_eigvals_p=0.0", 20, verbose=True
    ... )
    Saved data with 20 realizations to data/mc_eigvals_p=0.0_na=20.pkl
    """
    fname = working_path / f"{fname_base}_na={current_na}.pkl"
    with open(fname, "wb") as f:
        pk.dump(data, f)

    if metadata is not None:
        metadata_fname = working_path / f"{fname_base}_na={current_na}_metadata.pkl"
        with open(metadata_fname, "wb") as f:
            pk.dump(metadata, f)

    if verbose:
        print(f"Saved data with {current_na} realizations to {fname}")

    return fname


def cleanup_intermediate_files(
    working_path: Path,
    fname_base: str,
    final_na: int,
    verbose: bool = False
) -> None:
    """
    Remove all intermediate files with na < final_na.

    This prevents accumulation of partial result files after computation completes.

    Parameters
    ----------
    working_path : Path
        Directory containing files to clean.
    fname_base : str
        Base filename pattern (e.g., "mc_eigvals_p=0.0").
    final_na : int
        Final number of averages (files with na < this are removed).
    verbose : bool, optional
        Enable logging output.

    Returns
    -------
    None

    Examples
    --------
    >>> cleanup_intermediate_files(
    ...     Path("data/"), "mc_eigvals_p=0.0", 20, verbose=True
    ... )
    Removed intermediate file: mc_eigvals_p=0.0_na=10.pkl
    """
    glob_pattern = f"{fname_base}_na=*.pkl"
    for fpath in working_path.glob(glob_pattern):
        try:
            fname = fpath.stem
            na_part = fname.split('na=')[1]
            na_val = int(na_part)
            if na_val < final_na:
                fpath.unlink()
                # Also remove metadata if exists
                metadata_path = fpath.parent / f"{fpath.stem}_metadata.pkl"
                if metadata_path.exists():
                    metadata_path.unlink()
                if verbose:
                    print(f"Removed intermediate file: {fpath.name}")
        except (IndexError, ValueError):
            continue


def process_eigen_distribution(
    fname_base: str,
    initial_data_fn: Callable[[Any], np.ndarray],
    update_data_fn: Callable[[int, np.ndarray, np.ndarray, Any, Any], Any],
    save_data_fn: Callable[[Any, Any, Path], None],
    args: Any,
    get_working_path_fn: Callable[[Any], Path],
) -> None:
    """
    Manage computation and incremental saving of eigenvalue/eigenvector distributions.

    This generic function handles:
    - Checking for existing data (exact match or resume from partial)
    - Loading partial results
    - Batch processing with periodic saves
    - Automatic resume from interruption
    - Cleanup of intermediate files

    Parameters
    ----------
    fname_base : str
        Base filename for saving (without directory or extension).
        Example: "mc_dist_eigval_p=0.1" or "dist5_p=0.05_smallest"
    initial_data_fn : Callable[[Any], np.ndarray]
        Function to compute initial data: `initial_data_fn(args) -> np.ndarray`
    update_data_fn : Callable
        Function to update bins with signature:
        `update_data_fn(batch_size, bins, bin_centers, bin_counter, args) -> bin_counter`
    save_data_fn : Callable[[Any, Any, Path], None]
        Function to save data: `save_data_fn(args, data, filepath) -> None`
    args : Any
        Argument object with attributes:
        - `number_of_averages` : int - Total averages needed
        - `save_frequency` : int - Save checkpoint every N iterations
        - `verbose` : bool - Enable logging
        - `mode` : str - "eigval_dist" or "eigvec_dist"
        - `bins_count` : int - Number of histogram bins
        - `howmany` : int - Number of eigenvectors (for eigvec_dist)
    get_working_path_fn : Callable[[Any], Path]
        Function to get save directory: `get_working_path_fn(args) -> Path`

    Returns
    -------
    None
        Results saved to disk at `working_path / f"{fname_base}_na={navg_done}.pkl"`

    Examples
    --------
    >>> def my_get_path(args):
    ...     return Path("data/spectrum")
    >>> def my_initial_fn(args):
    ...     return compute_eigenvalues(args.L, args.p)
    >>> def my_update_fn(batch_size, bins, bin_centers, counter, args):
    ...     # Compute batch, update counter
    ...     return counter
    >>>
    >>> process_eigen_distribution(
    ...     "test_eigval_p=0.05",
    ...     my_initial_fn,
    ...     my_update_fn,
    ...     save_data,
    ...     args,
    ...     my_get_path
    ... )
    """
    if args.verbose:
        print(f"Starting process_eigen_distribution for {fname_base}")

    # Get working directory from callback
    working_path = get_working_path_fn(args)

    # Use new helper to find existing data
    existing_file, navg_done, should_skip = find_existing_data_by_na(
        working_path, fname_base, args.number_of_averages, args.verbose
    )

    if should_skip:
        return

    # Load existing data or compute initial data if no files exist
    if navg_done > 0 and existing_file is not None:
        bin_counter, initial_data = load_existing_data(
            existing_file, args.mode, args.verbose
        )
    else:
        initial_data = initial_data_fn(args)
        if initial_data is None:
            raise ValueError(
                "Initial data is empty. Please check the input parameters or data generation."
            )
        bin_counter = create_initial_bin_counter(args.mode, getattr(args, 'howmany', 0))

    # Create bins based on mode
    if args.mode == "eigvec_dist":
        bins, bin_centers = create_symmetric_log_bins(initial_data, args.bins_count)
    else:
        bins, bin_centers = linear_binning_hist(initial_data, args.bins_count)

    nAvgNeed = args.number_of_averages - navg_done
    total_periods = (
        nAvgNeed // args.save_frequency
        + bool(nAvgNeed % args.save_frequency)
    )

    # Process each period until reaching the required number of averages
    last_saved_na = navg_done  # Track last checkpoint to clean up old ones
    for current_period in range(total_periods):
        if args.verbose:
            print(f"Processing period {current_period + 1}/{total_periods}")
        batch_size = min(
            nAvgNeed - current_period * args.save_frequency,
            args.save_frequency,
        )
        bin_counter = update_data_fn(batch_size, bins, bin_centers, bin_counter, args)
        navg_done += batch_size

        # Use helper to save with proper na suffix
        save_with_na(bin_counter, working_path, fname_base, navg_done, verbose=args.verbose)

        # Clean up previous checkpoint (keep only latest) if not final
        if navg_done < args.number_of_averages and last_saved_na > 0 and last_saved_na < navg_done:
            prev_file = working_path / f"{fname_base}_na={last_saved_na}.pkl"
            prev_metadata = working_path / f"{fname_base}_na={last_saved_na}_metadata.pkl"
            if prev_file.exists():
                prev_file.unlink()
                if args.verbose:
                    print(f"Removed previous checkpoint: {prev_file.name}")
            if prev_metadata.exists():
                prev_metadata.unlink()

        last_saved_na = navg_done

    # Clean up last checkpoint if different from final
    if last_saved_na > 0 and last_saved_na < args.number_of_averages:
        prev_file = working_path / f"{fname_base}_na={last_saved_na}.pkl"
        prev_metadata = working_path / f"{fname_base}_na={last_saved_na}_metadata.pkl"
        if prev_file.exists():
            prev_file.unlink()
            if args.verbose:
                print(f"Removed last checkpoint: {prev_file.name}")
        if prev_metadata.exists():
            prev_metadata.unlink()
