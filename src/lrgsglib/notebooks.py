from lrgsglib.config.funcs import move_to_rootf

from .shared import *
from .core import *
from .plotlib import *
from .utils.ipy import *

# Engine-agnostic graph factories (preferred interface in notebooks).
from .graphs import (
    SignedGraph,
    Lattice2D,
    Lattice3D,
    LatticeND,
    ErdosRenyi,
    StochasticBlockModel,
    HierarchicalModular,
    MultiplicativeCascade,
    VicsekGraph,
    DiracCombGraph,
    DiracBrushGraph,
    SierpinskiGraph,
    BarabasiAlbert,
    WattsStrogatz,
    FullyConnected,
    kRegularGraph,
    BipartiteGraph,
    RandomGeometric,
    ConfigurationModel,
    LFRBenchmark,
    HolmeKim,
    DualBarabasiAlbert,
    ExtendedBarabasiAlbert,
    DGMgraph,
)

# LRG / spectral utilities most commonly needed alongside SignedGraph methods.
from .utils.lrg import (
    get_graph_lspectrum,
    compute_entropy_observables_from_eigenvalues,
    specific_heat_tau_window,
    lapl_dists,
    extract_ultrametric_matrix,
    MakeLinkageMatrix,
    compute_normalized_linkage,
    compute_optimal_threshold,
    circular_layout_by_cluster,
    log_dendrogram,
    dendrogram_leaf_node_colors,
    compute_signed_diffusion_distance,
    compute_eigenmode_sign_distance,
    agmon_geodesic_distance,
)
from .graphs.nx.funcs.spectral import signed_laplacian_matrix

# Most-used library defaults (symbolic, reusable).
from .config.const import (
    LRSG_ENTROPY_STEP,
    DEFAULT_ENTROPY_LEXPONENT,
    DEFAULT_ENTROPY_HEXPONENT,
    L2D_SIDE1,
    L2D_GEO_SQR,
    L2D_P_C_DICT,
)

LINKAGE_METHOD = 'average'   # UPGMA — classical LRG choice
CMAP_CLUSTERS  = 'tab20'

move_to_rootf()
plt.style.use('ipynb/nb_plotsheet.mplstyle')
warnings.filterwarnings('ignore', category=scipy.sparse.SparseEfficiencyWarning)