DEFAULT_HOWMANY_EIGS = 0
DEFAULT_L2DSSPECT_MODE = 'eigvec_dist'

# MCG_SlaplSpect defaults
DEFAULT_MCGSSPECT_MODE = 'eigval_dist'

# Backend for eigendecomposition (scipy/numpy/cupy)
DEFAULT_SLAPLSPECT_BACKEND = 'scipy'  # Conservative default; cluster will use 'cupy'

# Sparse eigendecomposition strategy
DEFAULT_KEEP_SPARSE = None  # Auto-select based on N and sparsity