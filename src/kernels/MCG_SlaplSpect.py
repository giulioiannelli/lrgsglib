from lrgsglib import *
from lrgsglib.utils.basic.probability import create_symmetric_log_bins, linear_binning_hist
from lrgsglib.config.funcs import bin_eigenvalues
import numpy as np
import random
import pickle as pk
import os
import glob
from pathlib import Path
from collections import Counter

# ============================================================================
# Helper functions for MultiplicativeCascadeGraph spectral analysis
# ============================================================================

def eigv_for_mc_graph(p1, p2, p3, p4, fraction, iterations,
                      stochastic=False, periodic=False, variant='exp_clocks',
                      mode='all', pflip=0.0, seed=None):
    """
    Compute eigenvalues of signed Laplacian for MultiplicativeCascadeGraph.

    Parameters
    ----------
    p1, p2, p3, p4 : float
        Cascade seed probabilities
    fraction : float
        Node sampling fraction
    iterations : int
        Cascade iterations
    stochastic : bool
        Stochastic cascade mode
    periodic : bool
        Periodic boundary conditions
    variant : str
        'standard' or 'exp_clocks'
    mode : str
        'full' or 'some_k' where k is number of eigenvalues
    pflip : float
        Probability of flipping edge signs
    seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    eigenvalues : np.ndarray
        Signed Laplacian eigenvalues
    """
    from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    # Generate graph
    G = MultiplicativeCascadeGraph(
        p1=p1, p2=p2, p3=p3, p4=p4,
        fraction=fraction,
        iterations=iterations,
        stochastic=stochastic,
        periodic=periodic,
        variant=variant,
        pflip=pflip
    )

    # Compute eigenvalues based on mode
    if mode == 'full':
        G.compute_laplacian_spectrum_weigV()
        eigvals = G.eigv  # Full spectrum
    elif mode.startswith('some'):
        k = int(mode.split('_')[1])
        G.compute_k_eigvV(k=k)
        eigvals = G.eigv  # Partial spectrum (first k eigenvalues)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return eigvals


def eigV_for_mc_graph_ptch(p1, p2, p3, p4, fraction, iterations,
                            stochastic=False, periodic=False, variant='exp_clocks',
                            mode='smallest', howmany=5, pflip=0.0,
                            backend='scipy', seed=None):
    """
    Compute eigenvectors of signed Laplacian for MultiplicativeCascadeGraph.

    Parameters
    ----------
    p1, p2, p3, p4 : float
        Cascade seed probabilities
    fraction : float
        Node sampling fraction
    iterations : int
        Cascade iterations
    stochastic : bool
        Stochastic cascade mode
    periodic : bool
        Periodic boundary conditions
    variant : str
        'standard' or 'exp_clocks'
    mode : str
        'smallest', 'largest', or 'middle'
    howmany : int
        Number of eigenvectors to compute
    pflip : float
        Probability of flipping edge signs
    backend : str
        Computational backend: 'numpy', 'scipy', or 'cupy'
    seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    eigenvectors : np.ndarray
        Shape (howmany, N) where N is number of nodes
    """
    from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    G = MultiplicativeCascadeGraph(
        p1=p1, p2=p2, p3=p3, p4=p4,
        fraction=fraction,
        iterations=iterations,
        stochastic=stochastic,
        periodic=periodic,
        variant=variant,
        pflip=pflip
    )

    # Compute eigenvectors with specified backend
    if mode == 'smallest':
        G.compute_k_eigvV(backend=backend, k=howmany)
        eigvecs = G.eigV[:howmany, :]  # First howmany eigenvectors
    elif mode == 'largest':
        # Compute all, then take last howmany
        G.compute_laplacian_spectrum_weigV()
        eigvecs = G.eigV[-howmany:, :]
    elif mode == 'middle':
        # Compute all, then take middle howmany
        G.compute_laplacian_spectrum_weigV()
        N = G.N
        mid = N // 2
        start = mid - howmany // 2
        eigvecs = G.eigV[start:start+howmany, :]
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return np.abs(eigvecs)  # Shape (howmany, N)


# ============================================================================
# Main processing function (mirrors L2D_SlaplSpect)
# ============================================================================

def perform_spectral_calculations(args):
    """
    Main function to parse arguments and run spectral analysis.

    Supports three modes:
    1. 'eigval_dist': Eigenvalue distribution over multiple realizations
    2. 'eigvec_dist': Eigenvector component distribution
    3. 'eigvals': Raw eigenvalue computation and storage
    """

    # Determine mode and process accordingly
    if args.mode.endswith("dist"):
        if args.mode == "eigvec_dist":
            fname_base = f"mc_dist{args.howmany}_{args.pflip:.3g}_{args.eigen_mode}"
            initial_fn = eigvec_initial_data
            update_fn = eigvec_update_data
        elif args.mode == "eigval_dist":
            fname_base = f"mc_dist_eigval_{args.pflip:.3g}"
            initial_fn = eigval_initial_data
            update_fn = eigval_update_data

        # Process distribution
        process_eigen_distribution(
            fname_base, initial_fn, update_fn, save_data, args
        )

    elif args.mode == "eigvals":
        fname_base = f"mc_eigvals_{args.pflip:.3g}"

        # Determine eigenvalue mode: full spectrum or partial
        # If howmany <= 0, compute full spectrum (needed for entropy calculations)
        # If howmany > 0, compute first howmany eigenvalues
        eig_mode = 'full' if args.howmany <= 0 else f'some_{args.howmany}'

        eigvlist = []
        for idx in range(args.number_of_averages):
            eigv = eigv_for_mc_graph(
                p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
                fraction=args.fraction,
                iterations=args.iterations,
                stochastic=args.stochastic,
                periodic=args.periodic,
                variant=args.variant,
                mode=eig_mode,
                pflip=args.pflip
            )
            eigvlist.append(eigv)

            # Periodic saving
            if idx % args.period == 0:
                path_fname_base = get_mc_spectrum_path(args)
                path_fname = path_fname_base / Path(f"{fname_base}_{idx}.pkl")

                # Remove previous checkpoint
                if idx > 0:
                    path_fname_prev = path_fname_base / Path(f"{fname_base}_{idx - args.period}.pkl")
                    if os.path.exists(path_fname_prev):
                        os.remove(path_fname_prev)

                with open(path_fname, "wb") as f:
                    pk.dump(eigvlist, f)

        # Save final result after loop completes
        path_fname_base = get_mc_spectrum_path(args)
        final_idx = args.number_of_averages - 1
        path_fname_final = path_fname_base / Path(f"{fname_base}_{final_idx}.pkl")

        # Remove last checkpoint if different from final
        if final_idx % args.period != 0:
            last_checkpoint = (final_idx // args.period) * args.period
            path_fname_prev = path_fname_base / Path(f"{fname_base}_{last_checkpoint}.pkl")
            if os.path.exists(path_fname_prev):
                os.remove(path_fname_prev)

        with open(path_fname_final, "wb") as f:
            pk.dump(eigvlist, f)

        if args.verbose:
            print(f"Saved final eigenvalue list with {len(eigvlist)} realizations to {path_fname_final}")


# ============================================================================
# Distribution processing (reuse pattern from L2D)
# ============================================================================

def process_eigen_distribution(fname_base, initial_data_fn, update_data_fn,
                                save_data_fn, args):
    """
    Manages computation and saving of eigenvalue/eigenvector distributions.
    """
    if args.verbose:
        print(f"Starting process_eigen_distribution for {fname_base}")

    # Check for existing files and determine the number of averages already done
    navg_done = 0
    working_path = get_mc_spectrum_path(args)
    search_str = f"{fname_base}_*.pkl"
    existing_files = sorted(glob.glob(str(working_path / Path(search_str))))

    if existing_files:
        navg_done = max(int(os.path.splitext(f.split('_')[-1])[0])
                        for f in existing_files)

    # Load existing data or compute initial data if no files exist
    if navg_done > 0:
        path_fname = working_path / Path(f"{fname_base}_{navg_done}.pkl")
        with open(path_fname, "rb") as f:
            bin_counter = pk.load(f)
        if args.mode == "eigvec_dist":
            initial_data = [val for counter in bin_counter for val in counter.elements()]
        else:
            initial_data = list(bin_counter.elements())
        if args.verbose:
            print(f"Loaded existing data from {path_fname}")
    else:
        path_fname = working_path / Path(f"{fname_base}_{navg_done}.pkl")
        initial_data = initial_data_fn(args)
        if initial_data is None:
            raise ValueError("Initial data is empty. Please check the input parameters.")
        if args.mode == "eigvec_dist":
            bin_counter = [Counter() for _ in range(args.howmany)]
        else:
            bin_counter = Counter()

    # Create bins based on mode
    if args.mode == "eigvec_dist":
        bins, bin_centers = create_symmetric_log_bins(initial_data, args.bins_count)
    else:
        bins, bin_centers = linear_binning_hist(initial_data, args.bins_count)

    nAvgNeed = args.number_of_averages - navg_done
    total_periods = (nAvgNeed // args.period) + bool(nAvgNeed % args.period)

    # Process each period
    for current_period in range(total_periods):
        if args.verbose:
            print(f"Processing period {current_period + 1}/{total_periods}")
        batch_size = min(nAvgNeed - current_period * args.period, args.period)
        bin_counter = update_data_fn(batch_size, bins, bin_centers,
                                     bin_counter, args)
        save_data_fn(args, bin_counter, path_fname)
        navg_done += batch_size
        new_fname = working_path / Path(f"{fname_base}_{navg_done}.pkl")
        os.rename(path_fname, new_fname)
        path_fname = new_fname

    # Rename the final file
    fname_final = working_path / Path(f"{fname_base}_{navg_done}.pkl")
    os.rename(path_fname, fname_final)
    if args.verbose:
        print(f"Renamed final file to {fname_final}")


# ============================================================================
# Data generation functions
# ============================================================================

def eigvec_initial_data(args):
    """Compute initial eigenvector data for MC graph."""
    if args.verbose:
        print("Computing initial eigenvector data for MC graph...")

    result = eigV_for_mc_graph_ptch(
        p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
        fraction=args.fraction,
        iterations=args.iterations,
        stochastic=args.stochastic,
        periodic=args.periodic,
        variant=args.variant,
        mode=args.eigen_mode,
        howmany=args.howmany,
        pflip=args.pflip,
        backend=args.eigen_mode
    )

    if args.verbose:
        print(f"Computed eigenvector data of size {result.shape}")
    return result


def eigval_initial_data(args):
    """Compute initial eigenvalue data for MC graph."""
    if args.verbose:
        print("Computing initial eigenvalue data for MC graph...")

    result = eigv_for_mc_graph(
        p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
        fraction=args.fraction,
        iterations=args.iterations,
        stochastic=args.stochastic,
        periodic=args.periodic,
        variant=args.variant,
        mode='full',
        pflip=args.p
    )

    if args.verbose:
        print(f"Computed eigenvalue data of size {result.shape}")
    return result


def eigvec_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """Update bin counters for eigenvector values (batch processing)."""
    if args.verbose:
        print(f"Updating eigenvector data for batch size {batch_size}...")
    eig_values = [[] for _ in range(args.howmany)]

    for _ in range(batch_size):
        eigV = eigV_for_mc_graph_ptch(
            p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
            fraction=args.fraction,
            iterations=args.iterations,
            stochastic=args.stochastic,
            periodic=args.periodic,
            variant=args.variant,
            mode=args.eigen_mode,
            howmany=args.howmany,
            pflip=args.pflip,
            backend=args.eigen_mode
        )
        for i in range(args.howmany):
            eig_values[i].append(eigV[i])

    eig_values = [np.concatenate(i) for i in eig_values]

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
    """Update bin counters for eigenvalue values (batch processing)."""
    if args.verbose:
        print(f"Updating eigenvalue data for batch size {batch_size}...")

    eig_values = [eigv_for_mc_graph(
        p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
        fraction=args.fraction,
        iterations=args.iterations,
        stochastic=args.stochastic,
        periodic=args.periodic,
        variant=args.variant,
        mode='full',
        pflip=args.p
    ) for _ in range(batch_size)]

    eig_values = np.concatenate(eig_values)

    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = linear_binning_hist(eig_values, args.bins_count)

    bin_counter.update(bin_eigenvalues(eig_values, bins, bin_centers))

    if args.verbose:
        print("Updated eigenvalue data.")
    return bin_counter


def save_data(args, data, fname):
    """Save computed data to file."""
    if args.verbose:
        print(f"Saving data to {fname}...")
    with open(fname, "wb") as f:
        pk.dump(data, f)
    if args.verbose:
        print(f"Data saved to {fname}")


# ============================================================================
# Path utilities
# ============================================================================

def get_mc_spectrum_path(args):
    """
    Get the path for saving MC spectrum data.

    Returns
    -------
    path : Path
        Directory path for spectrum data storage
    """
    from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

    # Create a dummy graph to get path structure
    G = MultiplicativeCascadeGraph(
        p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
        fraction=args.fraction,
        iterations=args.iterations,
        only_const_mode=True  # Don't generate full graph
    )

    return G.path_spect
