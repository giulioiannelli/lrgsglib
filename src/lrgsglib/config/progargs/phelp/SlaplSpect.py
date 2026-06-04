phelp_howmany_eigs = "Number of eigenvalues/vectors to compute (use 0 or -1 for full spectrum in eigvals mode)"
phelp_l2dsspect_mode = "Mode of operation, either 'eigvec_dist' for \
    eigenvector distribution or 'eigval_dist' for eigenvalue distribution \
    or 'eigvals' for eigenvalues"

phelp_mcgsspect_mode = """
    Computation mode for MultiplicativeCascadeGraph spectral analysis:
    - 'eigval_dist': Eigenvalue distribution over multiple realizations
    - 'eigvec_dist': Eigenvector component distribution
    - 'eigvals': Compute and save eigenvalues for averaging
    - 'entropy': Compute entropy observables (bypasses eigenvalues via expm_multiply)
"""

phelp_entropy_steps = "Number of logarithmic time points for entropy profile (τ ∈ [10^t1, 10^t2])"
phelp_entropy_t1 = "Lower time exponent for entropy calculation (τ_min = 10^t1)"
phelp_entropy_t2 = "Upper time exponent for entropy calculation (τ_max = 10^t2)"
phelp_num_samples = "Number of Rademacher probes for Hutchinson trace estimator (error ∝ 1/√num_samples)"
phelp_entropy_seed = "Random seed for entropy probe vectors (None = non-deterministic)"
phelp_entropy_norm = "Entropy normalization: 'standard' (raw) or 'complement' (1 - S)"
phelp_specific_heat_scale = "Specific heat scaling: 'logN' (remove size dependence) or 'none'"

phelp_slaplspect_backend = """Backend for Laplacian eigendecomposition:
    - 'scipy': Optimized CPU (divide-and-conquer, ~2x faster for N>1000)
    - 'numpy': Standard CPU (compatibility)
    - 'cupy': GPU acceleration (requires CUDA)
    Default: 'scipy' (recommended for most cases)"""

phelp_keep_sparse = """Sparse eigendecomposition strategy for full spectrum:
    - None (auto): Use sparse for N>5000 with sparsity<0.3, else dense
    - True: Force sparse methods (computes N-2 eigenpairs, saves memory)
    - False: Force dense methods (computes all N eigenpairs)
    Default: None (auto-select)"""

phelp_laplacian_type = """Which signed Laplacian to diagonalize (Kunegis 2010):
    - 'signed': combinatorial L = D_s - A (symmetric)            [default]
    - 'sym':    symmetric normalized L_sym = I - D_s^-1/2 A D_s^-1/2
    - 'rw':     random-walk L_rw = I - D_s^-1 A (non-symmetric; isospectral
                to 'sym' but with the true random-walk eigenvectors)
    D_s is the signed (absolute) degree matrix. Default: 'signed'."""