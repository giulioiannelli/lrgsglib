"""
Signed Laplacian spectral analysis for Multiplicative Cascade Graphs.

This module provides spectral computation for MultiplicativeCascadeGraph:
- Eigenvalue distributions over stochastic realizations
- Eigenvector component distributions
- Raw eigenvalue storage

Uses generic framework from SlaplSpect.py with MCG-specific
graph construction functions.
"""

from lrgsglib import *
from lrgsglib.utils.basic.probability import create_symmetric_log_bins, linear_binning_hist
from lrgsglib.config.funcs import bin_eigenvalues
import numpy as np
import random
import pickle as pk
from pathlib import Path
from collections import Counter

from .SlaplSpect import (
    process_eigen_distribution,
    save_data,
    build_eigval_fname_base,
    build_eigvec_fname_base,
    find_existing_data_by_na,
    save_with_na,
    cleanup_intermediate_files,
)

# ============================================================================
# Helper functions for MultiplicativeCascadeGraph spectral analysis
# ============================================================================

def select_optimal_backend(N, E, requested_backend='cupy', verbose=False):
    """
    Select optimal backend and sparsity mode based on graph size.

    Prioritizes Dense GPU when possible, falls back to faster CPU method when needed.

    Parameters
    ----------
    N : int
        Number of nodes
    E : int
        Number of edges
    requested_backend : str
        User-requested backend ('cupy', 'scipy', 'numpy')
    verbose : bool
        Print selection reasoning

    Returns
    -------
    tuple[str, bool, str]
        (backend, keep_sparse, backend_suffix) where backend_suffix describes the choice
    """
    # Calculate memory requirements
    dense_gb = (N ** 2 * 8) / (1024 ** 3)
    sparse_mb = (E * 16) / (1024 ** 2)

    # A100 VRAM limit (conservative: 60GB out of 80GB to account for workspace)
    # Dense eigenvalue decomposition needs extra workspace (~2-3x matrix size)
    GPU_VRAM_LIMIT_GB = 60.0

    # CPU RAM limit (conservative, assumes ~200GB available on cluster nodes)
    CPU_RAM_LIMIT_GB = 150.0

    # Threshold where sparse becomes faster than dense on CPU
    # From benchmarks: sparse is slower for N < 50k, faster for N > 100k
    SPARSE_CROSSOVER_N = 75000

    # Strategy 1: Dense GPU if requested and fits
    if requested_backend == 'cupy' and dense_gb < GPU_VRAM_LIMIT_GB:
        if verbose:
            print(f"Backend auto-select: Dense GPU (cupy) - {dense_gb:.1f}GB < {GPU_VRAM_LIMIT_GB}GB limit")
        return ('cupy', False, f'denseGPU')

    # Strategy 2: Dense CPU if fits in RAM and not too large
    if dense_gb < CPU_RAM_LIMIT_GB and N < SPARSE_CROSSOVER_N:
        if verbose:
            print(f"Backend auto-select: Dense CPU (scipy) - {dense_gb:.1f}GB RAM, N={N:,} < {SPARSE_CROSSOVER_N:,}")
        return ('scipy', False, f'denseCPU')

    # Strategy 3: Sparse CPU for very large graphs
    # Gets N-2 eigenvalues (missing 2 is negligible for N > 10k)
    if verbose:
        print(f"Backend auto-select: Sparse CPU (scipy) - N={N:,} too large for dense")
        print(f"  Dense would be: {dense_gb:.1f}GB")
        print(f"  Sparse uses: {sparse_mb:.1f}MB")
        print(f"  Will compute N-2 = {N-2:,} eigenvalues (missing 2/{N} = {2/N*100:.4f}%)")
    return ('scipy', True, f'sparseCPU_N-2')


def eigv_for_mc_graph(p1, p2, p3, p4, fraction, iterations,
                      stochastic=False, periodic=False, variant='exp_clocks',
                      mode='all', pflip=0.0, backend='scipy', seed=None,
                      out_suffix='', keep_sparse=None, verbose=False):
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
    backend : str
        Computational backend: 'numpy', 'scipy', or 'cupy'
        If 'cupy', will auto-select optimal backend based on graph size
    seed : int, optional
        Random seed for reproducibility
    out_suffix : str
        Optional suffix for output directory
    keep_sparse : bool, optional
        Sparse strategy: None (auto), True (force sparse), False (force dense)
    verbose : bool
        If True, print backend selection messages

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
        pflip=pflip,
        out_suffix=out_suffix
    )

    # Auto-select optimal backend if cupy requested and keep_sparse not specified
    actual_backend = backend
    actual_keep_sparse = keep_sparse
    backend_suffix = ''

    if backend == 'cupy' and keep_sparse is None:
        # Use intelligent backend selection
        N = G.N
        E = G.gr['G'].number_of_edges()
        actual_backend, actual_keep_sparse, backend_suffix = select_optimal_backend(
            N, E, requested_backend=backend, verbose=verbose
        )
        # Append backend info to out_suffix if not already present
        if backend_suffix and backend_suffix not in out_suffix:
            G.out_suffix = f"{out_suffix}_{backend_suffix}" if out_suffix else backend_suffix

    # Compute eigenvalues based on mode
    if mode == 'full':
        # Use eigenvalues-only computation (faster, no eigenvectors)
        # Try GPU, fall back to CPU if OOM
        try:
            G.compute_laplacian_spectrum(backend=actual_backend, keep_sparse=actual_keep_sparse)
            eigvals = G.eigv  # Full spectrum (or N-2 if sparse)
        except Exception as e:
            # Check if it's a GPU OOM error
            if 'cupy' in str(type(e).__module__) and 'OutOfMemory' in str(type(e).__name__):
                print(f"GPU Out of Memory! Falling back to Dense CPU")
                print(f"  Error: {e}")
                # Fall back to dense CPU
                actual_backend = 'scipy'
                actual_keep_sparse = False
                backend_suffix = 'denseCPU_fallback'
                G.out_suffix = f"{out_suffix}_{backend_suffix}" if out_suffix else backend_suffix
                G.compute_laplacian_spectrum(backend=actual_backend, keep_sparse=actual_keep_sparse)
                eigvals = G.eigv
            else:
                # Re-raise if not OOM
                raise
    elif mode.startswith('some'):
        k = int(mode.split('_')[1])
        G.compute_k_eigvV(k=k, backend=actual_backend)
        eigvals = G.eigv  # Partial spectrum (first k eigenvalues)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return eigvals


def eigV_for_mc_graph_ptch(p1, p2, p3, p4, fraction, iterations,
                            stochastic=False, periodic=False, variant='exp_clocks',
                            mode='smallest', howmany=5, pflip=0.0,
                            backend='scipy', seed=None, out_suffix=''):
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
    out_suffix : str
        Optional suffix for output directory

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
        pflip=pflip,
        out_suffix=out_suffix
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


def get_mc_spectrum_path(args):
    """
    Get the path for saving MC spectrum data.

    Parameters
    ----------
    args : Any
        Argument object with MCG parameters.

    Returns
    -------
    path : Path
        Directory path for spectrum data storage
    """
    from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

    # Get out_suffix from args if present
    out_suffix = getattr(args, 'out_suffix', '')

    # Create a dummy graph to get path structure
    G = MultiplicativeCascadeGraph(
        p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
        fraction=args.fraction,
        iterations=args.iterations,
        stochastic=getattr(args, 'stochastic', False),
        out_suffix=out_suffix,
        only_const_mode=True  # Don't generate full graph
    )

    return G.path_spect


# ============================================================================
# Main processing function (uses generic framework)
# ============================================================================

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
        Argument object with MCG parameters and processing mode.

    Returns
    -------
    None
        Results saved to disk.
    """

    # Determine mode and process accordingly
    if args.mode.endswith("dist"):
        if args.mode == "eigvec_dist":
            fname_base = build_eigvec_fname_base(
                "mc_", args.howmany, args.pflip, args.eigen_mode
            )
            initial_fn = eigvec_initial_data
            update_fn = eigvec_update_data
        elif args.mode == "eigval_dist":
            fname_base = build_eigval_fname_base("mc_", args.pflip)
            initial_fn = eigval_initial_data
            update_fn = eigval_update_data

        # Process distribution using generic framework
        process_eigen_distribution(
            fname_base, initial_fn, update_fn, save_data, args, get_mc_spectrum_path
        )

    elif args.mode == "eigvals":
        # Raw eigenvalue storage mode
        # Use "p=" prefix for pflip to distinguish from matrix probabilities
        from lrgsglib.config.const import DEFAULT_P_FSTR_FMT
        from lrgsglib.config.progargs.defs.SlaplSpect import DEFAULT_HOWMANY_EIGS
        pflip_str = f"p={args.pflip:{DEFAULT_P_FSTR_FMT}}"

        # Determine eigenvalue mode: full spectrum or partial
        # In eigvals mode, default to full spectrum (needed for entropy calculations)
        # User can limit with --howmany > 1 if they want partial spectrum
        # If howmany is default (1) or <= 0, compute full spectrum
        # If howmany > 1, compute first howmany eigenvalues
        if args.howmany <= DEFAULT_HOWMANY_EIGS:
            # Default behavior for eigvals mode: compute full spectrum
            eig_mode = 'full'
        else:
            # User explicitly requested partial spectrum
            eig_mode = f'some_{args.howmany}'

        # Build filename base (without na suffix)
        fname_base = f"mc_eigvals_{pflip_str}"

        # Get out_suffix
        out_suffix = getattr(args, 'out_suffix', '')

        # Get working path
        working_path = get_mc_spectrum_path(args)

        # Use helper to check for existing data
        existing_file, start_idx, should_skip = find_existing_data_by_na(
            working_path, fname_base, args.number_of_averages, args.verbose
        )

        if should_skip:
            return

        # Load existing data or start fresh
        eigvlist = []
        metadata_list = []

        if start_idx > 0 and existing_file is not None:
            try:
                with open(existing_file, "rb") as f:
                    eigvlist = pk.load(f)
                metadata_file = existing_file.parent / f"{existing_file.stem}_metadata.pkl"
                if metadata_file.exists():
                    with open(metadata_file, "rb") as f:
                        metadata_list = pk.load(f)
            except Exception as e:
                if args.verbose:
                    print(f"Warning: Could not load existing file {existing_file}: {e}")
                    print("Starting from scratch")
                eigvlist = []
                metadata_list = []
                start_idx = 0

        # Compute remaining realizations
        last_saved_na = start_idx  # Track last checkpoint to clean up old ones
        for idx in range(start_idx, args.number_of_averages):
            eigv = eigv_for_mc_graph(
                p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
                fraction=args.fraction,
                iterations=args.iterations,
                stochastic=args.stochastic,
                periodic=args.periodic,
                variant=args.variant,
                mode=eig_mode,
                pflip=args.pflip,
                backend=args.backend,
                out_suffix=out_suffix,
                keep_sparse=getattr(args, 'keep_sparse', None),
                verbose=(idx == start_idx and args.verbose)
            )
            eigvlist.append(eigv)

            # Store metadata for this realization
            metadata_list.append({
                'realization_idx': idx,
                'p1': args.p1,
                'p2': args.p2,
                'p3': args.p3,
                'p4': args.p4,
                'fraction': args.fraction,
                'iterations': args.iterations,
                'stochastic': args.stochastic,
                'periodic': args.periodic,
                'variant': args.variant,
                'pflip': args.pflip,
                'num_nodes': len(eigv),
                'backend': args.backend,
                'out_suffix': out_suffix,
                'eig_mode': eig_mode
            })

            # Periodic saving: save with current count to allow recovery
            if (idx + 1) % args.save_frequency == 0 and (idx + 1) < args.number_of_averages:
                current_count = idx + 1
                save_with_na(eigvlist, working_path, fname_base, current_count,
                           metadata=metadata_list, verbose=args.verbose)

                # Clean up previous checkpoint (keep only latest)
                if last_saved_na > 0 and last_saved_na < current_count:
                    prev_file = working_path / f"{fname_base}_na={last_saved_na}.pkl"
                    prev_metadata = working_path / f"{fname_base}_na={last_saved_na}_metadata.pkl"
                    if prev_file.exists():
                        prev_file.unlink()
                        if args.verbose:
                            print(f"Removed previous checkpoint: {prev_file.name}")
                    if prev_metadata.exists():
                        prev_metadata.unlink()

                last_saved_na = current_count

        # Save final result using helper
        save_with_na(eigvlist, working_path, fname_base, args.number_of_averages,
                   metadata=metadata_list, verbose=args.verbose)

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


# ============================================================================
# Data generation functions (MCG-specific implementations)
# ============================================================================

def eigvec_initial_data(args):
    """
    Compute initial eigenvector data for MC graph.

    Parameters
    ----------
    args : Any
        Argument object with MCG parameters.

    Returns
    -------
    np.ndarray
        Absolute values of eigenvector components.
    """
    if args.verbose:
        print("Computing initial eigenvector data for MC graph...")

    out_suffix = getattr(args, 'out_suffix', '')

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
        backend=args.backend,
        out_suffix=out_suffix
    )

    if args.verbose:
        print(f"Computed eigenvector data of size {result.shape}")
    return result


def eigval_initial_data(args):
    """
    Compute initial eigenvalue data for MC graph.

    Parameters
    ----------
    args : Any
        Argument object with MCG parameters.

    Returns
    -------
    np.ndarray
        Eigenvalues of signed Laplacian.
    """
    if args.verbose:
        print("Computing initial eigenvalue data for MC graph...")

    out_suffix = getattr(args, 'out_suffix', '')

    result = eigv_for_mc_graph(
        p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
        fraction=args.fraction,
        iterations=args.iterations,
        stochastic=args.stochastic,
        periodic=args.periodic,
        variant=args.variant,
        mode='full',
        pflip=args.pflip,
        out_suffix=out_suffix,
        verbose=args.verbose
    )

    if args.verbose:
        print(f"Computed eigenvalue data of size {result.shape}")
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
        Argument object with MCG parameters.

    Returns
    -------
    list of Counter
        Updated bin counters.
    """
    if args.verbose:
        print(f"Updating eigenvector data for batch size {batch_size}...")
    eig_values = [[] for _ in range(args.howmany)]

    out_suffix = getattr(args, 'out_suffix', '')

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
            backend=args.backend,
            out_suffix=out_suffix
        )
        for i in range(args.howmany):
            eig_values[i].append(eigV[i])

    eig_values = [np.concatenate(i) for i in eig_values]

    # Dynamically adjust bins if values fall outside current range
    min_val = min(arr.min() for arr in eig_values)
    max_val = max(arr.max() for arr in eig_values)
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
        Argument object with MCG parameters.

    Returns
    -------
    Counter
        Updated bin counter.
    """
    if args.verbose:
        print(f"Updating eigenvalue data for batch size {batch_size}...")

    out_suffix = getattr(args, 'out_suffix', '')

    eig_values = [eigv_for_mc_graph(
        p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
        fraction=args.fraction,
        iterations=args.iterations,
        stochastic=args.stochastic,
        periodic=args.periodic,
        variant=args.variant,
        mode='full',
        pflip=args.pflip,
        out_suffix=out_suffix,
        verbose=(i == 0 and args.verbose)
    ) for i in range(batch_size)]

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
