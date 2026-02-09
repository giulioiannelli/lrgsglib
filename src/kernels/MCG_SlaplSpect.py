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
from parsers.shared import get_graph_engine, resolve_graph_class

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

    # Detect GPU memory if cupy is available
    GPU_VRAM_LIMIT_GB = 70.0  # Default for A100-80GB cluster
    if requested_backend == 'cupy':
        try:
            import cupy as cp
            gpu_mem_bytes = cp.cuda.Device(0).mem_info[1]  # Total memory
            gpu_mem_gb = gpu_mem_bytes / (1024 ** 3)
            # Use 87.5% of available GPU memory (leave ~10-12.5% for overhead)
            GPU_VRAM_LIMIT_GB = gpu_mem_gb * 0.875
        except:
            pass  # Fall back to default if cupy not available

    # CPU RAM limit for dense eigendecomposition
    # High-memory nodes (Nix/Orion) have 768-1536 GB, so we can handle larger matrices
    CPU_RAM_LIMIT_GB = 500.0

    # NOTE: This function is called for FULL SPECTRUM computation
    # Sparse ARPACK with k≈N allocates N×N workspace → always prefer dense!
    # Sparse is ONLY useful for partial spectrum (k << N) - not handled here

    # Strategy 1: Dense GPU if requested and fits
    if requested_backend == 'cupy' and dense_gb < GPU_VRAM_LIMIT_GB:
        if verbose:
            print(f"Backend auto-select: Dense GPU (cupy) - {dense_gb:.1f}GB < {GPU_VRAM_LIMIT_GB:.1f}GB limit", flush=True)
        return ('cupy', False, f'denseGPU')

    # Strategy 2: Dense CPU if fits in RAM
    # For full spectrum, dense is ALWAYS better than sparse (no N×N workspace issue)
    if dense_gb < CPU_RAM_LIMIT_GB:
        if verbose:
            print(f"Backend auto-select: Dense CPU (scipy) - {dense_gb:.1f}GB < {CPU_RAM_LIMIT_GB:.1f}GB limit", flush=True)
            print(f"  Full spectrum for N={N:,} nodes", flush=True)
        return ('scipy', False, f'denseCPU')

    # Strategy 3: Graph too large for full spectrum
    # User must either use partial spectrum (--howmany) or reduce graph size
    if verbose:
        print(f"Backend auto-select: Graph too large for full spectrum!", flush=True)
        print(f"  N={N:,} nodes would require {dense_gb:.1f}GB (limit: {CPU_RAM_LIMIT_GB:.1f}GB)", flush=True)
        print(f"  Solution: Use --howmany to compute partial spectrum", flush=True)
    # Return sparse as placeholder, but it will likely fail
    return ('scipy', True, f'tooLarge_usePartialSpectrum')


def _make_mc_graph(p1, p2, p3, p4, fraction, iterations,
                   stochastic=False, periodic=False, variant='exp_clocks',
                   pflip=0.0, out_suffix='', seed=None, engine='nx'):
    """Create a MultiplicativeCascadeGraph instance for the given engine."""
    if engine == 'gt':
        Cls = resolve_graph_class("MultiplicativeCascade", "gt")
    else:
        from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as Cls
    return Cls(
        p1=p1, p2=p2, p3=p3, p4=p4,
        fraction=fraction,
        iterations=iterations,
        stochastic=stochastic,
        periodic=periodic,
        variant=variant,
        pflip=pflip,
        out_suffix=out_suffix,
        seed=seed,
    )


def eigv_for_mc_graph(p1, p2, p3, p4, fraction, iterations,
                      stochastic=False, periodic=False, variant='exp_clocks',
                      mode='all', pflip=0.0, backend='scipy', seed=None,
                      out_suffix='', keep_sparse=None, verbose=False,
                      engine='nx'):
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
    engine : str
        Graph engine: 'nx' (NetworkX) or 'gt' (graph-tool)

    Returns
    -------
    eigenvalues : np.ndarray
        Signed Laplacian eigenvalues
    """
    # Memory tracking for debugging OOM issues
    if verbose:
        try:
            import psutil
            import os
            process = psutil.Process(os.getpid())
            def get_mem_gb():
                return process.memory_info().rss / (1024**3)
            mem_start = get_mem_gb()
            print(f"[MEM] Start: {mem_start:.2f} GB", flush=True)
        except ImportError:
            verbose_mem = False
            print("[WARNING] psutil not available, memory tracking disabled", flush=True)
    else:
        verbose_mem = False

    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    # Generate graph
    if verbose:
        print(f"[STEP 1/5] Creating MultiplicativeCascadeGraph (p1={p1}, p2={p2}, p3={p3}, p4={p4}, fraction={fraction}, it={iterations}, engine={engine})...", flush=True)

    G = _make_mc_graph(
        p1=p1, p2=p2, p3=p3, p4=p4,
        fraction=fraction,
        iterations=iterations,
        stochastic=stochastic,
        periodic=periodic,
        variant=variant,
        pflip=pflip,
        out_suffix=out_suffix,
        seed=seed,
        engine=engine,
    )

    if verbose:
        try:
            mem_after_graph = get_mem_gb()
            print(f"[MEM] After graph creation: {mem_after_graph:.2f} GB (+{mem_after_graph - mem_start:.2f} GB)", flush=True)
            print(f"[INFO] Graph: N={G.N:,} nodes, E={G.Ne:,} edges", flush=True)
        except:
            pass

    # Auto-select optimal backend if cupy requested and keep_sparse not specified
    actual_backend = backend
    actual_keep_sparse = keep_sparse
    backend_suffix = ''

    if backend == 'cupy' and keep_sparse is None:
        # Use intelligent backend selection
        N = G.N
        E = G.Ne
        if verbose:
            print(f"[DEBUG] Graph has N={N:,} nodes, E={E:,} edges", flush=True)
        actual_backend, actual_keep_sparse, backend_suffix = select_optimal_backend(
            N, E, requested_backend=backend, verbose=verbose
        )
        # Append backend info to out_suffix if not already present
        if backend_suffix and backend_suffix not in out_suffix:
            G.out_suffix = f"{out_suffix}_{backend_suffix}" if out_suffix else backend_suffix

    # Compute eigenvalues based on mode
    if mode == 'full':
        # Use eigenvalues-only computation (faster, no eigenvectors)
        # CRITICAL: For full spectrum (k=N or k=N-2), MUST use dense method
        # Sparse ARPACK with k≈N allocates N×N workspace, defeating the purpose!
        if actual_keep_sparse is None:
            actual_keep_sparse = False  # Force dense for full spectrum
            if verbose:
                print(f"[INFO] Forcing keep_sparse=False for full spectrum (sparse with k≈N is inefficient)", flush=True)

        if verbose:
            print(f"[STEP 2/5] Computing Laplacian spectrum (backend={actual_backend}, keep_sparse={actual_keep_sparse})...", flush=True)
            matrix_gb = (G.N ** 2 * 8) / (1024**3)
            print(f"[INFO] Dense matrix size: {matrix_gb:.1f} GB", flush=True)

        try:
            if hasattr(G, 'compute_laplacian_spectrum'):
                G.compute_laplacian_spectrum(backend=actual_backend, keep_sparse=actual_keep_sparse, verbose=verbose)
            else:
                # GT graphs only have compute_laplacian_spectrum_weigV (no keep_sparse)
                G.compute_laplacian_spectrum_weigV(backend=actual_backend)
            eigvals = G.eigv  # Full spectrum

            if verbose:
                try:
                    mem_after_spectrum = get_mem_gb()
                    print(f"[MEM] After spectrum computation: {mem_after_spectrum:.2f} GB (+{mem_after_spectrum - mem_after_graph:.2f} GB)", flush=True)
                    print(f"[INFO] Computed {len(eigvals)} eigenvalues", flush=True)
                except:
                    pass
        except Exception as e:
            # Check if it's a GPU error that we can recover from
            error_str = str(e).lower()
            error_type = str(type(e))
            is_gpu_oom = 'cupy' in str(type(e).__module__) and 'OutOfMemory' in str(type(e).__name__)
            is_cusolver_limit = 'CUSOLVERError' in error_type or 'CUSOLVER_STATUS_INVALID_VALUE' in str(e)

            if is_gpu_oom or is_cusolver_limit:
                if is_gpu_oom:
                    print(f"GPU Out of Memory! Falling back to Dense CPU", flush=True)
                else:
                    print(f"cuSolver limit exceeded (N={G.N:,} too large for 32-bit API)! Falling back to Dense CPU", flush=True)
                print(f"  Error: {e}", flush=True)
                # Fall back to dense CPU
                actual_backend = 'scipy'
                actual_keep_sparse = False
                backend_suffix = 'denseCPU_fallback'
                G.out_suffix = f"{out_suffix}_{backend_suffix}" if out_suffix else backend_suffix
                if hasattr(G, 'compute_laplacian_spectrum'):
                    G.compute_laplacian_spectrum(backend=actual_backend, keep_sparse=actual_keep_sparse, verbose=verbose)
                else:
                    G.compute_laplacian_spectrum_weigV(backend=actual_backend)
                eigvals = G.eigv
            else:
                # Re-raise if not a recoverable GPU error
                raise
    elif mode.startswith('some'):
        k = int(mode.split('_')[1])
        G.compute_k_eigvV(k=k, backend=actual_backend)
        eigvals = G.eigv  # Partial spectrum (first k eigenvalues)
    else:
        raise ValueError(f"Unknown mode: {mode}", flush=True)

    return eigvals


def eigV_for_mc_graph_ptch(p1, p2, p3, p4, fraction, iterations,
                            stochastic=False, periodic=False, variant='exp_clocks',
                            mode='smallest', howmany=5, pflip=0.0,
                            backend='scipy', seed=None, out_suffix='',
                            engine='nx'):
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
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    G = _make_mc_graph(
        p1=p1, p2=p2, p3=p3, p4=p4,
        fraction=fraction,
        iterations=iterations,
        stochastic=stochastic,
        periodic=periodic,
        variant=variant,
        pflip=pflip,
        out_suffix=out_suffix,
        seed=seed,
        engine=engine,
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
        raise ValueError(f"Unknown mode: {mode}", flush=True)

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
    from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph

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
    out_suffix = getattr(args, 'out_suffix', '')
    engine = get_graph_engine(args)

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

        # Append out_suffix to fname_base if provided
        if out_suffix:
            fname_base = f"{fname_base}_{out_suffix}"

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
        # Include out_suffix in filename if provided
        if out_suffix:
            fname_base = f"mc_eigvals_{pflip_str}_{out_suffix}"
        else:
            fname_base = f"mc_eigvals_{pflip_str}"

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
                    print(f"Warning: Could not load existing file {existing_file}: {e}", flush=True)
                    print("Starting from scratch", flush=True)
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
                verbose=(idx == start_idx and args.verbose),
                engine=engine,
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
                            print(f"Removed previous checkpoint: {prev_file.name}", flush=True)
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
                    print(f"Removed last checkpoint: {prev_file.name}", flush=True)
            if prev_metadata.exists():
                prev_metadata.unlink()

    elif args.mode == "entropy":
        # Entropy mode: compute entropy observables from eigenvalues
        # Uses full spectrum decomposition (much faster than expm_multiply for tau ranges up to 10^7)

        from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues
        from lrgsglib.config.const import DEFAULT_P_FSTR_FMT

        pflip_str = f"p={args.pflip:{DEFAULT_P_FSTR_FMT}}"

        # Build filename base with out_suffix if provided
        if out_suffix:
            fname_base = f"mc_entropy_{pflip_str}_{out_suffix}"
        else:
            fname_base = f"mc_entropy_{pflip_str}"

        # Get entropy-specific parameters from args
        entropy_steps = getattr(args, 'entropy_steps', 600)
        entropy_t1 = getattr(args, 'entropy_t1', -2)
        entropy_t2 = getattr(args, 'entropy_t2', 5)
        entropy_seed = getattr(args, 'entropy_seed', None)
        entropy_norm = getattr(args, 'entropy_norm', 'complement')
        specific_heat_scale = getattr(args, 'specific_heat_scale', 'logN')
        spectral_backend = getattr(args, 'spectral_backend', 'scipy')  # scipy, numpy, or cupy

        # Get working path
        working_path = get_mc_spectrum_path(args)

        # Check for existing data
        existing_file, start_idx, should_skip = find_existing_data_by_na(
            working_path, fname_base, args.number_of_averages, args.verbose
        )

        if should_skip:
            return

        # Load existing data or start fresh
        entropy_list = []
        tau_scale = None  # Will be set from first computation

        if start_idx > 0 and existing_file is not None:
            try:
                with open(existing_file, "rb") as f:
                    data = pk.load(f)
                    if isinstance(data, tuple) and len(data) == 2:
                        entropy_list, tau_scale = data
                    else:
                        # Old format fallback
                        entropy_list = data if isinstance(data, list) else []
                if args.verbose:
                    print(f"Loaded {len(entropy_list)} existing entropy profiles", flush=True)
            except Exception as e:
                if args.verbose:
                    print(f"Warning: Could not load {existing_file}: {e}", flush=True)
                    print("Starting from scratch", flush=True)
                entropy_list = []
                tau_scale = None
                start_idx = 0

        # Compute remaining realizations
        last_saved_na = start_idx

        for idx in range(start_idx, args.number_of_averages):
            if entropy_seed is not None:
                # Use deterministic seed sequence
                realization_seed = entropy_seed + idx
            else:
                realization_seed = None

            if realization_seed is not None:
                np.random.seed(realization_seed)
                random.seed(realization_seed)

            G = _make_mc_graph(
                p1=args.p1, p2=args.p2, p3=args.p3, p4=args.p4,
                fraction=args.fraction,
                iterations=args.iterations,
                stochastic=args.stochastic,
                periodic=args.periodic,
                variant=args.variant,
                pflip=args.pflip,
                out_suffix=out_suffix,
                seed=realization_seed,
                engine=engine,
            )

            # Compute Laplacian spectrum
            if hasattr(G, 'compute_laplacian_spectrum'):
                G.compute_laplacian_spectrum(backend=spectral_backend)
            else:
                G.compute_laplacian_spectrum_weigV(backend=spectral_backend)

            if G.eigv is None:
                raise ValueError(f"Eigenvalue computation failed for realization {idx}", flush=True)

            # Compute entropy observables from eigenvalues
            entropy, specific_heat, variance, time_grid = compute_entropy_observables_from_eigenvalues(
                eigenvalues=G.eigv,
                num_nodes=G.N,
                steps=entropy_steps,
                t1=entropy_t1,
                t2=entropy_t2,
                entropy_norm=entropy_norm,
                specific_heat_scale=specific_heat_scale
            )

            if idx == start_idx and args.verbose:
                print(f"  Computed entropy for N={G.N} nodes", flush=True)
                print(f"  Eigenvalue range: [{np.min(G.eigv):.6f}, {np.max(G.eigv):.6f}]", flush=True)
                print(f"  Tau range: [{time_grid[0]:.2e}, {time_grid[-1]:.2e}]", flush=True)

            # Store tau_scale from first realization (same for all)
            if tau_scale is None:
                tau_scale = time_grid

            # Accumulate results
            entropy_list.append({
                'entropy': entropy,
                'specific_heat': specific_heat,
                'variance': variance,
                'num_nodes': G.N,
                'realization_idx': idx
            })

            # Periodic saving
            if (idx + 1) % args.save_frequency == 0 and (idx + 1) < args.number_of_averages:
                current_count = idx + 1
                save_with_na((entropy_list, tau_scale), working_path, fname_base, current_count,
                           metadata=None, verbose=args.verbose)

                # Clean up previous checkpoint
                if last_saved_na > 0 and last_saved_na < current_count:
                    prev_file = working_path / f"{fname_base}_na={last_saved_na}.pkl"
                    if prev_file.exists():
                        prev_file.unlink()
                        if args.verbose:
                            print(f"Removed previous checkpoint: {prev_file.name}", flush=True)

                last_saved_na = current_count

        # Save final result
        save_with_na((entropy_list, tau_scale), working_path, fname_base, args.number_of_averages,
                   metadata=None, verbose=args.verbose)

        # Clean up last checkpoint
        if last_saved_na > 0 and last_saved_na < args.number_of_averages:
            prev_file = working_path / f"{fname_base}_na={last_saved_na}.pkl"
            if prev_file.exists():
                prev_file.unlink()
                if args.verbose:
                    print(f"Removed last checkpoint: {prev_file.name}", flush=True)

        if args.verbose:
            print(f"\nEntropy mode completed:", flush=True)
            print(f"  - Computed {len(entropy_list)} realizations", flush=True)
            print(f"  - Time points: {len(tau_scale)}", flush=True)
            print(f"  - Saved to: {fname_base}_na={args.number_of_averages}.pkl", flush=True)


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
        print("Computing initial eigenvector data for MC graph...", flush=True)

    out_suffix = getattr(args, 'out_suffix', '')

    engine = get_graph_engine(args)
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
        out_suffix=out_suffix,
        engine=engine,
    )

    if args.verbose:
        print(f"Computed eigenvector data of size {result.shape}", flush=True)
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
        print("Computing initial eigenvalue data for MC graph...", flush=True)

    out_suffix = getattr(args, 'out_suffix', '')

    engine = get_graph_engine(args)
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
        verbose=args.verbose,
        engine=engine,
    )

    if args.verbose:
        print(f"Computed eigenvalue data of size {result.shape}", flush=True)
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
        print(f"Updating eigenvector data for batch size {batch_size}...", flush=True)
    eig_values = [[] for _ in range(args.howmany)]

    out_suffix = getattr(args, 'out_suffix', '')

    engine = get_graph_engine(args)
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
            out_suffix=out_suffix,
            engine=engine,
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
        print("Updated eigenvector data.", flush=True)
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
        print(f"Updating eigenvalue data for batch size {batch_size}...", flush=True)

    out_suffix = getattr(args, 'out_suffix', '')

    engine = get_graph_engine(args)
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
        verbose=(i == 0 and args.verbose),
        engine=engine,
    ) for i in range(batch_size)]

    eig_values = np.concatenate(eig_values)

    # Dynamically adjust bins if values fall outside current range
    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = linear_binning_hist(eig_values, args.bins_count)

    bin_counter.update(bin_eigenvalues(eig_values, bins, bin_centers))

    if args.verbose:
        print("Updated eigenvalue data.", flush=True)
    return bin_counter
