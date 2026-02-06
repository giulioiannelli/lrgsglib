"""
Signed Laplacian spectral analysis for 2D lattices.

This module provides spectral computation for Lattice2D graphs:
- Eigenvalue distributions over disorder realizations
- Eigenvector component distributions
- Raw eigenvalue storage for entropy analysis

Uses generic framework from SlaplSpect.py with Lattice2D-specific
graph construction functions.
"""

from lrgsglib import *
from .L2D import *
from .SlaplSpect import (
    process_eigen_distribution,
    save_data,
    build_eigval_fname_base,
    build_eigvec_fname_base,
)


def get_lattice2d_path(args):
    """
    Get spectrum save path for Lattice2D.

    Parameters
    ----------
    args : Any
        Argument object with L, p, geo, and workDir attributes.

    Returns
    -------
    Path
        Directory path for spectrum data storage.
    """
    return Lattice2D(
        side1=args.L, pflip=args.p, geo=args.geo, path_data=args.workDir
    ).path_spect


def perform_spectral_calculations(args):
    """
    Main function to parse arguments and run spectral analysis.

    Supports three modes:
    1. 'eigval_dist': Eigenvalue distribution over multiple realizations
    2. 'eigvec_dist': Eigenvector component distribution
    3. 'eigvals': Raw eigenvalue computation and storage

    Parameters
    ----------
    args : Any
        Argument object with mode, L, p, geo, number_of_averages, etc.

    Returns
    -------
    None
        Results saved to disk.
    """
    # Determine the filename base and processing functions based on the mode
    if args.mode.endswith("dist"):
        if args.mode == "eigvec_dist":
            fname_base = build_eigvec_fname_base(
                "", args.howmany, args.p, args.eigen_mode
            )
            initial_fn = eigvec_initial_data
            update_fn = eigvec_update_data
        elif args.mode == "eigval_dist":
            fname_base = build_eigval_fname_base(
                "", args.p, args.cell_type
            )
            initial_fn = eigval_initial_data
            update_fn = eigval_update_data

        # Process the eigen distribution using generic framework
        process_eigen_distribution(
            fname_base,
            initial_fn,
            update_fn,
            save_data,
            args,
            get_lattice2d_path
        )

    elif args.mode == "eigvals":
        # Raw eigenvalue storage mode
        fname_base = f"eigvals_{args.p:.3g}_{args.cell_type}"

        eigvlist = []
        for idx in range(args.number_of_averages):
            eigv = eigv_for_lattice2D(
                side=args.L,
                pflip=args.p,
                geo=args.geo,
                mode='_'.join(["some", str(args.howmany)]),
                backend=args.backend,
                keep_sparse=getattr(args, 'keep_sparse', None)
            )
            eigvlist.append(eigv)

            # Periodic saving with checkpoint management
            if idx % args.save_frequency == 0:
                path_fname_base = get_lattice2d_path(args)
                path_fname = path_fname_base / Path(f"{fname_base}_{idx}.pkl")

                # Remove previous checkpoint
                if idx > 0:
                    path_fname_prev = path_fname_base / Path(
                        f"{fname_base}_{idx - args.save_frequency}.pkl"
                    )
                    if os.path.exists(path_fname_prev):
                        os.remove(path_fname_prev)

                with open(path_fname, "wb") as f:
                    pk.dump(eigvlist, f)


def eigvec_initial_data(args):
    """
    Compute initial eigenvector data for Lattice2D.

    Parameters
    ----------
    args : Any
        Argument object with L, p, geo, eigen_mode, and howmany attributes.

    Returns
    -------
    np.ndarray
        Absolute values of eigenvector components.
    """
    if args.verbose:
        print("Computing initial eigenvector data...")
    result = np.abs(eigV_for_lattice2D_ptch(
        side=args.L,
        pflip=args.p,
        geo=args.geo,
        mode=args.eigen_mode,
        howmany=args.howmany
    ))
    if args.verbose:
        print(f"Computed initial eigenvector data of size {result.shape}")
    return result


def eigval_initial_data(args):
    """
    Compute initial eigenvalue data for Lattice2D.

    Parameters
    ----------
    args : Any
        Argument object with L, p, and geo attributes.

    Returns
    -------
    np.ndarray
        Eigenvalues of signed Laplacian.
    """
    if args.verbose:
        print("Computing initial eigenvalue data...")
    result = eigv_for_lattice2D(
        side=args.L,
        pflip=args.p,
        geo=args.geo
    )
    if args.verbose:
        print(f"Computed initial eigenvalue data of size {result.shape}")
    return result


def eigvec_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """
    Update bin counters for eigenvector values (batch processing).

    Parameters
    ----------
    batch_size : int
        Number of graph realizations to process in this batch.
    bins : np.ndarray
        Bin edges for histogram.
    bin_centers : np.ndarray
        Bin center values.
    bin_counter : list of Counter
        Current bin counters (one per eigenvector).
    args : Any
        Argument object with graph parameters.

    Returns
    -------
    list of Counter
        Updated bin counters.
    """
    if args.verbose:
        print(f"Updating eigenvector data for batch size {batch_size}...")

    eig_values = [[] for _ in range(args.howmany)]
    for _ in range(batch_size):
        eigV = eigV_for_lattice2D_ptch(
            side=args.L,
            pflip=args.p,
            geo=args.geo,
            mode=args.eigen_mode,
            howmany=args.howmany
        )
        for i in range(args.howmany):
            eig_values[i].append(eigV[i])

    eig_values = [np.concatenate(i) for i in eig_values]

    # Dynamically adjust bins if values fall outside current range
    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = create_symmetric_log_bins(eig_values, args.bins_count)

    for i in range(args.howmany):
        bin_counter[i].update(bin_eigenvalues(eig_values[i], bins, bin_centers))

    if args.verbose:
        print("Updated eigenvector data.")
    return bin_counter


def eigval_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """
    Update bin counters for eigenvalue values (batch processing).

    Parameters
    ----------
    batch_size : int
        Number of graph realizations to process in this batch.
    bins : np.ndarray
        Bin edges for histogram.
    bin_centers : np.ndarray
        Bin center values.
    bin_counter : Counter
        Current bin counter.
    args : Any
        Argument object with graph parameters.

    Returns
    -------
    Counter
        Updated bin counter.
    """
    if args.verbose:
        print(f"Updating eigenvalue data for batch size {batch_size}...")

    eig_values = [
        eigv_for_lattice2D(side=args.L, pflip=args.p, geo=args.geo)
        for _ in range(batch_size)
    ]
    eig_values = np.concatenate(eig_values)

    # Dynamically adjust bins if values fall outside current range
    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = linear_binning_hist(eig_values, args.bins_count)

    bin_counter.update(bin_eigenvalues(eig_values, bins, bin_centers))

    if args.verbose:
        print("Updated eigenvalue data.")
    return bin_counter
