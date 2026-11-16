#
DEFAULT_L2D_TRANSCLUSTER_MODE = "ordParam"
DEFAULT_L2D_TRANSCLUSTER_FLOAT_TYPE = "float64"
DEFAULT_L2D_TRANSCLUSTER_SAVE_FREQUENCY = 0
DEFAULT_L2D_TRANSCLUSTER_OUT_SUFFIX = ""
DEFAULT_L2D_TRANSCLUSTER_NAVG = 1000
# Default lists for serialiser batch submission L2D
DEFAULT_L2D_TRANSCLUSTER_SRUN_L_LIST = (10, 20, 30, 40)
DEFAULT_L2D_TRANSCLUSTER_SRUN_P_LIST = (0.1, 0.15, 0.2, 0.25)
DEFAULT_L2D_TRANSCLUSTER_SRUN_MU_LIST = None
DEFAULT_L2D_TRANSCLUSTER_SRUN_SIGMA_LIST = None
DEFAULT_L3D_TRANSCLUSTER_MODE = "ordParam"
DEFAULT_L3D_TRANSCLUSTER_FLOAT_TYPE = "float64"
DEFAULT_L3D_TRANSCLUSTER_SAVE_FREQUENCY = 0
DEFAULT_L3D_TRANSCLUSTER_OUT_SUFFIX = ""
DEFAULT_L3D_TRANSCLUSTER_NAVG = 1000
DEFAULT_L3D_TRANSCLUSTER_PDIL = 0.0
DEFAULT_L3D_TRANSCLUSTER_MU = 0.0
DEFAULT_L3D_TRANSCLUSTER_SIGMA = 0.0
DEFAULT_L3D_TRANSCLUSTER_EDGE_WEIGHT = "flip"
# Default lists for serialiser batch submission L3D
DEFAULT_L3D_TRANSCLUSTER_SRUN_L_LIST = (10, 20, 30, 40)
DEFAULT_L3D_TRANSCLUSTER_SRUN_P_LIST = (0.1, 0.15, 0.2, 0.25)
DEFAULT_L3D_TRANSCLUSTER_SRUN_MU_LIST = None
DEFAULT_L3D_TRANSCLUSTER_SRUN_SIGMA_LIST = None
# Default lists for serialiser batch submission SCS
DEFAULT_SCS_TRANSCLUSTER_SRUN_N_LIST = (8, 10, 12, 14)
DEFAULT_SCS_TRANSCLUSTER_SRUN_GAMMA_LIST = (0.1, 0.2, 0.3, 0.4)
DEFAULT_SCS_TRANSCLUSTER_SRUN_J0_LIST = (0.0, 0.5, 1.0)
DEFAULT_SCS_TRANSCLUSTER_SRUN_J_LIST = (1.0, 2.0, 3.0)
DEFAULT_SCS_TRANSCLUSTER_SRUN_G_LIST = (0.0, 0.1, 0.2)
# TransCluster v2 core (kernels/TransClusterCore.py): graph-agnostic nested
# disorder -> signed Laplacian ground state -> sign clusters -> record rows.
# Seeds: SeedSequence(root, spawn_key=(SCHEMA, *graph_key, rep, STREAM, SUB)).
TC_SEED_SCHEMA = 1  # first spawn-key word; bump if the seed tree changes
TC_SEED_SPAN = 2**32 - 1  # engine seeds live in [1, TC_SEED_SPAN]; 0 = entropy
TC_KEY_WORD_LIMIT = 2**32  # every spawn-key word must stay below this
TC_STREAM_TOPO = 0  # graph generator seed (per rep and attempt)
TC_STREAM_SIGNS = 1  # S family: u_e ~ U(0, 1) per edge
TC_STREAM_XNODES = 2  # X family: w_i ~ U(0, 1) per node
TC_STREAM_GAUSS = 3  # G family: z_e ~ N(0, 1) per edge
TC_STREAM_V0 = 4  # Lanczos start vector ~ N(0, 1)^N (never the ones vector)
TC_STREAM_SUB = 0  # reserved sub-stream word (always 0 in v1)
TC_FAMILY_SIGNS = "S"  # -c rand: p = M^-/M negative edges
TC_FAMILY_XNODES = "X"  # -c randXERR: q = fraction of X nodes, star edges SET
TC_FAMILY_GAUSS = "G"  # -ew normal: J_e = mu + sigma * z_e
TC_FAMILIES = (TC_FAMILY_SIGNS, TC_FAMILY_XNODES, TC_FAMILY_GAUSS)
TC_ENSEMBLE_CANONICAL = "canonical"  # exact count m = round-half-up(p * M)
TC_ENSEMBLE_BERNOULLI = "bernoulli"  # independent u_e < p
TC_ENSEMBLES = (TC_ENSEMBLE_CANONICAL, TC_ENSEMBLE_BERNOULLI)
DEFAULT_TC_ENSEMBLE = TC_ENSEMBLE_CANONICAL
DEFAULT_TC_MU = 1.0  # G family mean coupling (mu = 0 makes sigma irrelevant)
TC_DENSE_MAX_N = 256  # dense eigh at or below this N, Lanczos above
# Lanczos basis. At N = 1e6, band edge, it needs 2.2-2.8x fewer matvecs than
# 20 for the same signs (+155 MB); 80 saves only 1.2-1.7x more at N = 1e5.
TC_EIG_NCV = 40
# Matvec cap of one solve attempt, a deterministic cost bound (no clock).
# A unit makes at most 3 attempts (TransClusterCore.solve_attempts). At
# N = 1e6 a nev = 1 solve needs about 4-6k matvecs, a nev = 2 solve 12k+.
TC_EIG_MAX_MV = 20_000
TC_EIG_RETRY_NCV_FACTOR = 2  # the nev = 1 retry uses a larger basis
DEFAULT_TC_EIG_TOL = 0.0  # 0 = machine precision
# Largest --eig-tol the program accepts. At N = 2000, tol = 1e-3 changed P_inf
# in 21 of 48 units, 1e-4 in 1 and 1e-6 in none; tol >= 1 checks nothing.
TC_EIG_TOL_MAX = 1e-6
# A returned eigenvector must have unit norm within this tolerance, or the
# attempt failed: ARPACK can return x = 0 without raising.
TC_EIG_NORM_TOL = 1e-8
# Every residual ||L x - lambda x|| of a Lanczos attempt must be at most this
# times max(eps, tol) ||L||, or the attempt failed. Converged pairs measured
# 1-9 (66 solves, N = 2e4-1e5); pairs ARPACK accepts on couplings of scale
# 1e-20 to 1e-26, where its own test is absolute, 1e7 to 3e13.
TC_EIG_RESID_FACTOR = 1000.0
DEFAULT_TC_NEV = 1  # eigenpairs per solve (gap reps solve TC_GAP_NEV)
TC_GAP_NEV = 2  # gap reps add lambda_1, its eigenvector and energy
DEFAULT_TC_GAP_EVERY = 10  # reps with rep % n == 0 are gap reps (0 = off)
# |x_i| < c * eps * ||L|| is below the resolution limit
# (iannelli2025topological.pdf, SI Sec. IV.A); ||L|| <= 2 max_i sum_j |J_ij|.
TC_NUNREL_EPS_FACTOR = 1.0
TC_GERSHGORIN_FACTOR = 2.0
TC_FINGERPRINT_BYTES = 8  # blake2b digest size of the topology fingerprint
TC_STATUS_OK = 0
# No attempt converged within TC_EIG_MAX_MV matvecs with its residuals
# inside TC_EIG_RESID_FACTOR.
TC_STATUS_NOCONV = 1
TC_STATUS_NONFINITE = 2  # NaN, inf or a vector of norm != 1 came back
TC_STATUS_ERROR = 3  # the solver raised (an ARPACK error code, LinAlgError)
TC_NO_COUNT = -1  # cnt column when no count is prescribed (bernoulli, G)
# CLI vocabulary of the v2 programs: -c picks S or X, -ew normal picks G.
TC_CELL_FAMILY = {"rand": TC_FAMILY_SIGNS, "randXERR": TC_FAMILY_XNODES}
DEFAULT_TC_CELL = "rand"
TC_EDGE_WEIGHT_FLIP = "flip"  # +-1 couplings, family from -c
TC_EDGE_WEIGHT_NORMAL = "normal"  # gaussian couplings (family G)
TC_EDGE_WEIGHTS = (TC_EDGE_WEIGHT_FLIP, TC_EDGE_WEIGHT_NORMAL)
DEFAULT_TC_EDGE_WEIGHT = TC_EDGE_WEIGHT_FLIP
TC_LANG = "np"  # the lang= run-dir token: pure numpy/scipy pipeline
# Allowed control values per family: fractions p, q in [0, 1]; sigma >= 0.
TC_CTRL_RANGE = {
    TC_FAMILY_SIGNS: (0.0, 1.0),
    TC_FAMILY_XNODES: (0.0, 1.0),
    TC_FAMILY_GAUSS: (0.0, float("inf")),
}
# G coupling scale max(|mu|, sigma). L scales with it, so the physics does not
# depend on it, but the numerics do: ARPACK's convergence test turns absolute
# when ||L|| nears eps^(2/3) ~ 4e-11 (P_inf changed at sigma <= 1e-25, mu = 0),
# and residual norms overflow from about 1e170, degree sums from about 3e307.
TC_G_SCALE_MIN = 1e-8
TC_G_SCALE_MAX = 1e100
TC_TAG_FORBIDDEN = "_=/"  # -os tag: no run-dirname separators, no paths
TC_RUNDIR_NAME_MAX = (
    255  # bytes (UTF-8); NAME_MAX of the usual Linux file systems
)
# Tiles (kernels/TransClusterTiles.py): a tile is reps [r0, r1) x control
# indices [c0, c1); tile ids decode rep-major.
TC_TILE_SCHEMA = 1  # bump if the tile npz layout changes
TC_GRID_HASH_BYTES = 4  # blake2b digest of the grid: 8 hex chars in names
DEFAULT_TC_REP_START = 0
DEFAULT_TC_CTRL_CHUNK = 0  # control values per tile (0 = all)
DEFAULT_TC_REP_CHUNK = 10  # reps per tile (0 = all)
TC_ALL_TILES = -1  # tile id that runs every tile in this process
DEFAULT_TC_TILE_ID = TC_ALL_TILES
TC_CKPT_SECONDS = 300.0  # checkpoint at the first unit boundary after this
TC_STOP_SIGNALS = ("SIGUSR1", "SIGTERM")  # checkpoint, then exit
TC_TILES_DIR = "tiles"
TC_TILE_EXT = ".npz"
TC_PART_EXT = ".part.npz"
TC_BAD_EXT = ".bad"  # suffix of a tile file set aside (unreadable or foreign)
TC_BAD_ID_BYTES = 4  # random id in a set-aside name: 8 hex characters
TC_SLURM_ENV_VARS = (
    "SLURM_JOB_ID",
    "SLURM_ARRAY_JOB_ID",
    "SLURM_ARRAY_TASK_ID",
    "SLURM_JOB_NODELIST",
)
TC_PROVENANCE_MODULES = ("numpy", "scipy", "networkx", "graph_tool")
TC_GIT_TIMEOUT = 10.0  # seconds for each provenance git call
TC_EXIT_OK = 0
TC_EXIT_CFG_MISMATCH = 3  # cfg.json or a tile belongs to another run
TC_EXIT_PREEMPTED = 4  # signalled: partial written, resubmit the tile
TC_EXIT_FAILED = (
    1  # serializer: sbatch or its log failed, or nothing to collect
)
# Serializer (RRG_TransCluster_Serializer.py): SLURM arrays of tiles.
TC_SPLIT_CTRL_MIN_N = 100_000  # from this N: TC_SPLIT_CTRL_CHUNK per tile
TC_SPLIT_CTRL_CHUNK = 1  # one control value per tile
DEFAULT_TC_ARRAY_CONCURRENT = 200  # the %K throttle of each array (not a total)
DEFAULT_TC_ARRAY_PARTITION = ""  # blank = cluster default
DEFAULT_TC_ARRAY_TIME = "04:00:00"
DEFAULT_TC_MEM_MB = 0  # 0 = the estimate below
TC_SIGNAL_LEAD_SECONDS = 300  # stop signal this long before the time limit
# Memory estimate per task (MB), measured on a laptop on 2026-09-30:
# SAFETY * (BASE + N * (PER_NODE + ENGINE[-ge])). Recalibrate on the cluster.
TC_MEM_BASE_MB = 450.0
TC_MEM_MB_PER_NODE = 4e-4
TC_MEM_ENGINE_MB_PER_NODE = {"gt": 7.6e-4, "nx": 1.5e-3}
TC_MEM_SAFETY = 1.5
TC_SUBMISSIONS_DIR = "submissions"  # one <job id>.json per submitted array
# Statistics (kernels/TransClusterStats.py)
TC_AGG_SCHEMA = 1  # bump if the agg npz layout changes
TC_AGG_NAME = "agg.npz"
TC_LEGACY_TXT_NAME = "ordParam.txt"
# One multi-row file sorted by control value (legacy: one file per p).
# Columns: converged reps, control value (sigma for -ew normal), <n_neg>,
# <P>, <P^2>, sd(P, ddof=0).
TC_LEGACY_FMT = "%.7g"
TC_KURT_MIN_N = 10  # fewer converged reps: kurtosis and its chi SE are NaN
DEFAULT_TC_BOOT = 1000  # paired bootstrap resamples of the chi peak
DEFAULT_TC_BOOT_SEED = 0
TC_BOOT_SEED_LIMIT = 2**63  # --boot-seed lies in [0, this)
TC_BOOT_QUANTILES = (0.16, 0.84)  # 68 % interval of the chi peak
# Fewer paired reps: the chi-peak interval is NaN. n reps have only
# C(2n-1, n) distinct resamples (3 at n = 2), so the interval collapses.
TC_BOOT_MIN_REPS = 10
TC_PEAK_MIN_CTRLS = 3  # control values of the parabola through the chi max
TC_SRUN_ACTIONS = ("print", "exec", "status", "collect")  # serializer modes
DEFAULT_TC_SRUN_ACTION = "print"  # no side effects
