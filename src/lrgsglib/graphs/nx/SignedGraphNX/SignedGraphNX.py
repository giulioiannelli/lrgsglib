import copy
import os
import logging
import random
import struct
import time
import warnings

import cupy as cp
import networkx as nx
import numpy as np
import pickle as pk

from networkx import Graph
from typing import Any, Union, List, Dict, Tuple, Optional, Type, Callable
from numbers import Number
from matplotlib.pyplot import get_cmap
from numpy.typing import NDArray
from scipy.sparse import spdiags
from scipy.sparse import identity as scsp_identity
from scipy.sparse.linalg import eigsh as scsp_eigsh

from ....config.const import *
from ....config.funcs import build_p_fname
from ....config.errwar import NflipError, NoClustError, SignedGraphWarning
from ....utils.basic import is_in_range, join_non_empty,\
    flip_to_positive_majority_adapted, bin_sign, flip_to_positive_majority,\
    generate_random_id, normalize_array, dtype_numerical_precision
from ....utils.lrg import compute_ising_pairwise_energy
from ....utils.tools import NestedDict, ConditionalPartitioning
from ....utils.tools.ConditionalPartitioning import ConditionalPartitioningInput
from ._backend import Backend, BackendManager, ArrayBackend
#
logger = logging.getLogger(__name__)

class SignedGraphNX:
    """
    A class for representing and analyzing signed graphs with positive and 
    negative edge weights.
    
    SignedGraph extends NetworkX graph functionality to support signed networks,
    where edges can have positive or negative weights. The class provides 
    comprehensive tools for analyzing structural properties, computing spectral 
    characteristics of the signed Laplacian, clustering analysis, and studying 
    spin-glass-like dynamics on graphs.
    
    Attributes
    ----------
    gr : dict
        Dictionary of graph representations (alias for 
        graph_representation_dictionary).
    G : Graph
        The primary NetworkX graph object.
    on_g : str
        Identifier for the primary graph representation being used.
    pflip : float
        Probability of edge sign flipping (proportion of negative edges).
    seed : int
        Random seed for reproducibility.
    N : int
        Number of nodes in the graph (property).
    Ne : int
        Number of edges in the graph (property).
    Ne_n : int
        Number of negative edges in the graph (property).
    adj : scipy.sparse matrix
        Adjacency matrix (property).
    lap : scipy.sparse matrix
        Standard Laplacian matrix (property).
    slp : scipy.sparse matrix
        Signed Laplacian matrix (property).
    degm : scipy.sparse matrix
        Degree matrix (property).
    sdeg : scipy.sparse matrix
        Signed degree matrix (property).
    eigv : NDArray
        Eigenvalues of the signed Laplacian.
    eigV : NDArray
        Eigenvectors of the signed Laplacian (can be row or column major).
    eset : dict
        Dictionary mapping graph representations to sets of all edges.
    fleset : dict
        Dictionary mapping graph representations to sets of negative edges.
    lfeset : dict
        Dictionary mapping graph representations to sets of positive edges.
    nwContainer : Type[Any], optional
        Network container class to be defined by subclasses. If defined and 
        init_nw_dict is True, an instance will be created and stored in nwDict.
    
    Methods
    -------
    **Core Graph Operations**:

    - ``flip_sel_edges(links, name='weight', on_g=SG_REPR)``:
      Flip the sign of selected edges.
    - ``flip_random_fract_edges(pflip=None, on_g=SG_REPR)``:
      Randomly flip a fraction of edges based on pflip.
    - ``get_random_links(n=1, only_in='', on_g=SG_REPR)``:
      Get random edges from the graph (all, positive, or negative).

    **Spectral Analysis** (imported from ._spectral):

    - ``compute_laplacian_spectrum_weigV(backend='numpy')``:
      Compute full eigendecomposition of the signed Laplacian.
    - ``compute_k_eigvV(k=1, backend='scipy', which='SM')``:
      Compute the k smallest eigenvalues/eigenvectors.
    - ``get_eigV(which=0, binarize=False)``:
      Retrieve a specific eigenvector (optionally binarized).
    - ``get_eigV_binarized(which=0)``:
      Get binarized eigenvector (sign pattern).
    - ``get_signed_laplacian_embedding(k=2)``:
      Get k-dimensional spectral embedding.

    **Clustering and Community Detection**:

    - ``get_eigV_cluster_sizes(which=0, binarize=True)``:
      Compute cluster sizes from eigenvector-based partitioning.
    - ``get_cluster_distribution(which=0, binarize=True)``:
      Get distribution of cluster sizes.
    - ``make_clustersYN(k, val, on_g=SG_REPR)``:
      Create clusters based on node attribute values.
    - ``compute_pinf(which=0, val=None, on_g=SG_REPR)``:
      Compute infinite cluster probability (percolation order parameter).

    **Energy and Dynamics** (imported from ._dynamics):

    - ``compute_rbim_energy_eigV(which=0, use_gpu=False)``:
      Compute Random Bond Ising Model energy for binarized eigenvector.
    - ``compute_sksph_energy_eigV(which=0, use_gpu=False)``:
      Compute spherical Sherrington-Kirkpatrick model energy.
    - ``get_ferroAntiferro_regions(attr_str='s', on_g=SG_REPR)``:
      Identify ferromagnetic and antiferromagnetic regions.

    **Information Theory** (imported from ._infotheory):

    - ``compute_signed_laplacian_entropy(steps=600, t1=-2, t2=5)``:
      Compute entropy profile and related observables from spectrum.
    - ``get_entropy()``:
      Return cached normalized entropy profile.
    - ``get_specific_heat()``:
      Return entropy derivative (specific heat analog).
    
    **Import/Export** (via ._exports / ._loaders):

    - ``export_eigV_all(out_suffix='', ext=BIN)``:
      Export all eigenvectors to binary file.
    - ``load_eigV_all(fname='', load_mode=SG_LOAD_M)``:
      Load eigenvectors from binary file.
    - ``set_edgel_from_bin(file_path, mode='numpy', on_g=SG_REPR)``:
      Load edge list from binary file.

    **Graph Representations**:

    - ``upd_GraphRepr_All(on_g=SG_REPR, also_itself=True)``:
      Update all graph representations to maintain consistency.
    - ``upd_graph_matrices(format='csr', on_g=SG_REPR)``:
      Update adjacency, degree, and Laplacian matrices.
    
    Notes
    -----
    The signed Laplacian is defined as L_s = D_s - A, where D_s is the
    diagonal matrix of absolute degree values and A is the signed adjacency
    matrix. This formulation is particularly useful for studying frustration,
    balance theory, and spin-glass dynamics on networks [1]_ [2]_.

    Edge weights are expected to be ±1 or arbitrary real values. Negative
    edges represent antagonistic or frustrated interactions, while positive
    edges represent cooperative interactions.

    The class supports multiple graph representations that can be maintained
    simultaneously (e.g., different node labelings) through the `on_g`
    parameter.

    Spectral methods use either NumPy, SciPy, or CuPy (GPU-accelerated)
    backends depending on availability and problem size.
    
    Examples
    --------

    Create a signed graph from a NetworkX graph:

    >>> import networkx as nx  # doctest: +SKIP
    >>> G = nx.erdos_renyi_graph(100, 0.1)  # doctest: +SKIP
    >>> sg = SignedGraph(G, pflip=0.3, seed=42)  # doctest: +SKIP
    >>> print(f"Nodes: {sg.N}, Edges: {sg.Ne}")  # doctest: +SKIP

    Compute and analyze the signed Laplacian spectrum:

    >>> sg.compute_laplacian_spectrum_weigV()  # doctest: +SKIP
    >>> sg.compute_signed_laplacian_entropy()  # doctest: +SKIP
    >>> entropy = sg.get_entropy()  # doctest: +SKIP

    Analyze clustering from eigenvectors:

    >>> cluster_sizes = sg.get_eigV_cluster_sizes(which=0)  # doctest: +SKIP
    >>> sg.compute_pinf(which=0)  # doctest: +SKIP

    Compute spin-glass energy:

    >>> sg.compute_rbim_energy_eigV(which=0)  # doctest: +SKIP
    >>> energy = sg.get_rbim_energy_eigV(which=0)  # doctest: +SKIP
    
    See Also
    --------
    networkx.Graph : Base NetworkX graph class
    scipy.sparse.linalg.eigsh : Sparse eigenvalue solver
    
    References
    ----------
    .. [1] Cucuringu, M. (2016). "Sync-Rank: Robust Ranking, Constrained 
           Ranking and Rank Aggregation via Eigenvector and SDP 
           Synchronization." IEEE Transactions on Network Science and 
           Engineering.
    .. [2] Kirkpatrick, S., & Sherrington, D. (1978). "Infinite-ranged models 
           of spin-glasses." Physical Review B, 17(11), 4384.
    """
    sgpathn = "signed_graph"
    
    # Optional class-level attributes that can be overridden by subclasses
    nwContainer: Optional[Type[Any]] = None  # Network container class for subclasses
    syshape: Optional[int] = None  # System shape (set by topology subclasses)
    syshapePth: Optional[str] = None  # System shape path string
    #
    def __init__(
        self,
        G: Graph,
        pflip: float = SG_PFLIP,
        on_g: str = SG_REPR,
        seed: Optional[int] = None,
        path_data: Optional[Path] = None,
        path_plot: Optional[Path] = None,
        init_nw_dict: bool = SG_INIT_NW_DICT,
        init_weights_val: Union[float, Dict[str, float]] = SG_INIT_WVAL,
        export_mode: str = SG_EXPORT_M,
        make_dir_tree: bool = True,
        imported: bool = SG_IMPORT_ON,
        import_fname: str = '',
        import_mode: str = SG_LOAD_M,
        backend: Union[str, Backend] = Backend.NUMPY,
    ) -> None:
        """
        Initialize a SignedGraph instance.

        Parameters
        ----------
        G : Graph
            A NetworkX graph object to be used as the base graph.
        pflip : float, default SG_PFLIP
            Probability of flipping edges to negative weights (0.0 to 1.0).
        on_g : str, default SG_REPR
            Primary graph representation identifier.
        seed : int, optional
            Random seed for reproducibility. If None, generates from time/PID.
        path_data : Path, optional
            Base path for data storage. Defaults to global PATHDATA.
        path_plot : Path, optional
            Base path for plot storage. Defaults to global PATHPLOT.
        init_nw_dict : bool, default SG_INIT_NW_DICT
            Whether to initialize the network dictionary container.
        init_weights_val : Union[float, Dict[str, float]], default SG_INIT_WVAL
            Initial edge weight value(s) for new graphs.
        export_mode : str, default SG_EXPORT_M
            Export mode for graph data (currently unused in init).
        make_dir_tree : bool, default True
            Whether to create the directory structure on initialization.
        imported : bool, default SG_IMPORT_ON
            Whether the graph is being loaded from file.
        import_fname : str, default ''
            Filename to load graph from (if imported=True).
        import_mode : str, default SG_LOAD_M
            Mode for loading graph data.
        backend : Union[str, Backend], default Backend.NUMPY
            Computational backend for array operations. Can be 'numpy', 'scipy',
            'cupy', or a Backend enum value. CuPy enables GPU acceleration.

        Notes
        -----
        Initialization follows these steps:
        1. Initialize backend for numerical computations
        2. Initialize core data structures (edge sets, mappings)
        3. Validate pflip and set up randomness
        4. Ensure required attributes for subclass compatibility
        5. Initialize file paths and directory structure
        6. Set configuration parameters
        7. Load or assign the base graph
        8. Initialize graph representations
        9. Set up signed graph structure (unless only_const_mode)
        10. Optionally initialize network dictionary and clustering utility

        Subclasses (topology wrappers) may set attributes like `std_fname`,
        `only_const_mode`, and `nwContainer` before calling this constructor.
        """
        # Initialize backend for numerical computations
        self._init_backend(backend)

        # Initialize core data structures
        self.graph_representation_dictionary = {}
        self.map_node = {}
        self.map_edge = {}
        self.eset = {}
        self.fleset = {}
        self.lfeset = {}
        
        # Initialize matrix dictionaries (keyed by graph representation)
        self.adjacency_matrices = {}
        self.degree_matrices = {}
        self.signed_degree_matrices = {}
        self.laplacian_matrices = {}
        self.signed_laplacian_matrices = {}
        
        # Initialize clustering-related attributes
        self.clustersY: Optional[list[set]] = None
        self.clustersN: Optional[list[set]] = None
        self.numClustersY: int = 0
        self.numClustersN: int = 0
        self.biggestClSet: list[set] = []
        self.numClustersBig: int = 0
        self.gc: Optional[set] = None
        self.largest_cc: Optional[set] = None
        self.largest_cc_subgraph: Optional[Graph] = None
        
        # Validate pflip and initialize randomness
        self._verify_pflip(pflip)
        self.__init_randomness__(seed)
        
        # Ensure compatibility with topology subclasses
        if not hasattr(self, 'std_fname'):
            self.std_fname = 'sg'
        if not hasattr(self, 'only_const_mode'):
            self.only_const_mode = False
        
        # Initialize file system paths
        self.__init_paths__(path_data, path_plot, make_dir_tree)
        
        # Set configuration parameters
        self.on_g = on_g
        self.init_weights_val = init_weights_val
        self.load_g = imported
        self.init_nw_dict = init_nw_dict
        
        # Compose standard filename
        self.std_fname = join_non_empty('_', self.std_fname, self.peq_str)
        
        # Load or assign base graph
        self.G = (
            self.__load_graph__(import_fname, import_mode) 
            if self.load_g 
            else G
        )
        
        # Initialize graph representations
        self.__init_reprdict__()
        
        # Initialize signed graph structure (unless in const-only mode)
        if not self.only_const_mode:
            self.__init_sgraph__()
            
            # Initialize network dictionary if requested
            if self.init_nw_dict:
                if self.nwContainer is not None:
                    self.nwDict = self.nwContainer(self)
                else:
                    raise AttributeError(SG_ERRMSG_NW_DICT)
            
            # Initialize clustering utility
            self.__init_graph_clustering_utility__()

    def _init_backend(self, backend: Union[str, Backend]) -> None:
        """
        Initialize the computational backend for array operations.

        Parameters
        ----------
        backend : Union[str, Backend]
            Backend specification. Can be 'numpy', 'scipy', 'cupy', or a Backend enum value.

        Notes
        -----
        The backend is stored in two forms:
        - self._backend: The ArrayBackend instance for numerical operations
        - self._backend_name: String name of the backend ('numpy', 'scipy', or 'cupy')

        CuPy backend requires CUDA-capable GPU and CuPy installation.
        If CuPy is requested but unavailable, falls back to NumPy with warning.
        """
        # Store backend instance for numerical operations
        self._backend: ArrayBackend = BackendManager.get_backend(backend, fallback=True)
        self._backend_name: str = self._backend.name

    #
    @property
    def adj(self) -> Any:
        """Adjacency matrix of the signed graph for the current graph representation.

        Returns
        -------
        Any
            Sparse matrix (scipy.sparse) or None if not computed.
        """
        return self.adjacency_matrices.get(self.on_g)

    @property
    def degm(self) -> Any:
        """Degree matrix of the graph for the current graph representation.

        Returns
        -------
        Any
            Sparse diagonal matrix (scipy.sparse) or None if not computed.
        """
        return self.degree_matrices.get(self.on_g)

    @property
    def gr(self) -> Dict[str, Graph]:
        """Dictionary of graph representations.

        Returns
        -------
        Dict[str, Graph]
            Mapping from representation names to NetworkX graphs.
        """
        return self.graph_representation_dictionary

    @property
    def gcl(self) -> Any:
        """Graph clustering utility (nested dictionary structure).

        Returns
        -------
        Any
            Nested dictionary structure for clustering data.
        """
        return self.graph_clustering_utility

    @property
    def lap(self) -> Any:
        """Standard Laplacian matrix (D - A) for the current graph representation.

        Returns
        -------
        Any
            Sparse matrix (scipy.sparse) or None if not computed.
        """
        return self.laplacian_matrices.get(self.on_g)

    @property
    def sdeg(self) -> Any:
        """Signed degree matrix (absolute degree values) for the current graph representation.

        Returns
        -------
        Any
            Sparse diagonal matrix (scipy.sparse) or None if not computed.
        """
        return self.signed_degree_matrices.get(self.on_g)

    @property
    def slp(self) -> Any:
        """Signed Laplacian matrix (D_s - A) for the current graph representation.

        Returns
        -------
        Any
            Sparse matrix (scipy.sparse) or None if not computed.
        """
        return self.signed_laplacian_matrices.get(self.on_g)

    @property
    def N(self) -> int:
        """Number of nodes in the primary graph representation.

        Returns
        -------
        int
            Node count.
        """
        return self.gr[self.on_g].number_of_nodes()

    @property
    def Ne(self) -> int:
        """Number of edges in the primary graph representation.

        Returns
        -------
        int
            Edge count.
        """
        return self.gr[self.on_g].number_of_edges()

    @property
    def num_edges(self) -> int:
        """Number of edges (alias for Ne, GT compatibility)."""
        return self.gr[self.on_g].number_of_edges()

    @property
    def Ne_n(self) -> int:
        """Number of negative edges in the primary graph representation.

        Returns
        -------
        int
            Count of edges with negative weights.
        """
        return len(self.fleset[self.on_g])

    def count_negative_edges(self) -> int:
        """Count number of negative edges (GT compatibility).

        Returns
        -------
        int
            Number of edges with negative weight.
        """
        return len(self.fleset[self.on_g])

    def count_positive_edges(self) -> int:
        """Count number of positive edges (GT compatibility).

        Returns
        -------
        int
            Number of edges with positive weight.
        """
        return len(self.lfeset[self.on_g])

    def get_signed_adjacency(self) -> Any:
        """Get signed adjacency matrix (GT compatibility).

        Returns
        -------
        ndarray
            Signed adjacency matrix with +1/-1 edge weights.
        """
        return self.adj.toarray() if hasattr(self.adj, 'toarray') else self.adj

    def to_networkx(self) -> Graph:
        """Return the underlying NetworkX graph (GT compatibility).

        For NX-based graphs, this returns the primary graph representation.
        Edge weights are stored in the 'weight' attribute.

        Returns
        -------
        nx.Graph
            The NetworkX graph object.
        """
        return self.gr[self.on_g]

    @property
    def edges_with_data(self) -> Any:
        """Edges with data in the primary graph representation.

        Returns
        -------
        EdgeDataView
            NetworkX edge view with data attributes.
        """
        return self.gr[self.on_g].edges(data=True)

    def _invalidate_matrix_cache(self, on_g: Optional[str] = None) -> None:
        """
        Clear cached matrices after graph modification.

        This method should be called whenever the graph structure or edge weights
        are modified to ensure matrix properties return updated values.

        Parameters
        ----------
        on_g : str, optional
            Specific graph representation to invalidate. If None, invalidates
            all cached matrices for all representations.

        Notes
        -----
        The following matrices are cleared:
        - Adjacency matrices
        - Degree matrices
        - Signed degree matrices
        - Laplacian matrices
        - Signed Laplacian matrices

        Spectral properties (eigv, eigV) are NOT cleared as they may have been
        explicitly computed and should be managed separately.

        Examples
        --------

        >>> sg = SignedGraph(G=nx.karate_club_graph())  # doctest: +SKIP
        >>> sg.flip_random_fract_edges(pflip=0.2)  # doctest: +SKIP
        """
        if on_g is None:
            # Clear all cached matrices for all representations
            self.adjacency_matrices.clear()
            self.degree_matrices.clear()
            self.signed_degree_matrices.clear()
            self.laplacian_matrices.clear()
            self.signed_laplacian_matrices.clear()
        else:
            # Clear cached matrices for specific representation
            if on_g in self.adjacency_matrices:
                del self.adjacency_matrices[on_g]
            if on_g in self.degree_matrices:
                del self.degree_matrices[on_g]
            if on_g in self.signed_degree_matrices:
                del self.signed_degree_matrices[on_g]
            if on_g in self.laplacian_matrices:
                del self.laplacian_matrices[on_g]
            if on_g in self.signed_laplacian_matrices:
                del self.signed_laplacian_matrices[on_g]

    #
    def __ensure_required_attributes__(self) -> None:
        """
        Ensure all required attributes exist for backward compatibility
        with older pickled objects that may be missing some attributes.
        """
        # Initialize missing core attributes
        required_attrs = {
            'eset': {},
            'fleset': {},
            'lfeset': {},
            'map_node': {},
            'map_edge': {},
            'adjacency_matrices': {},
            'degree_matrices': {},
            'signed_degree_matrices': {},
            'laplacian_matrices': {},
            'signed_laplacian_matrices': {},
            # Will be initialized later if needed
            'graph_clustering_utility': None
        }
        
        for attr_name, default_value in required_attrs.items():
            if not hasattr(self, attr_name):
                setattr(self, attr_name, default_value)
    #
    def __init_randomness__(self, seed: Optional[int] = None) -> None:
        """
        Initialize random number generators with a reproducible seed.
        
        Sets the random seed for Python's `random`, NumPy, and CuPy (if 
        available) to ensure reproducibility across all random operations in 
        the signed graph. If no seed is provided, generates a pseudo-random 
        seed from the current time and process ID.
        
        Parameters
        ----------
        seed : int, optional
            Random seed value. If None, generates a seed from:
            `(time_us + process_id + object_id) mod (2^32 - 1)` to ensure 
            uniqueness across concurrent runs.
            
        Notes
        -----
        A random string identifier (`rand_str`) is also generated for 
        creating unique temporary filenames or identifiers.
        
        If CuPy is unavailable or fails to initialize, a warning is issued 
        but execution continues without GPU-accelerated random operations.
        """
        self.seed = seed or (
            (int(time.time() * 1_000_000) + os.getpid() + id(self)) % (2**32 - 1)
        )
        random.seed(self.seed)
        try:
            cp.random.seed(self.seed)
        except Exception:
            pass
        np.random.seed(self.seed)
        #
        self.rand_str = generate_random_id()
    #
    def __init_dirs__(self, exist_ok: bool = True):
        """
        Create the directory structure for storing signed graph data.
        
        Creates all directories specified in the `subpath_list` attribute, 
        which typically includes subdirectories for eigenvalues, eigenvectors, 
        Ising dynamics data, cluster information, and other computed 
        properties.
        
        Parameters
        ----------
        exist_ok : bool, default True
            If True, do not raise an error if directories already exist. 
            If False, raise FileExistsError for existing directories.
            
        Notes
        -----
        The directory structure is automatically populated during path 
        initialization in `__init_paths__()`. This method is called 
        internally when `make_dir_tree=True` during initialization.
        
        All directories are created relative to `path_sgdata`, which is 
        derived from the base data path and system shape parameters.
        """
        for _ in self.subpath_list: 
            os.makedirs(_, exist_ok=exist_ok)
    #
    def __init_graph_clustering_utility__(self):
        """
        Initialize the graph clustering utility as a nested dictionary.
        
        Creates a `NestedDict` structure to store clustering results from 
        eigenvector-based partitioning. The nested dictionary organizes 
        clustering information hierarchically by:
        - Attribute key (e.g., 'eigV0', 'eigV1')
        - Partitioning condition (e.g., '+1', '-1')
        - Graph representation (e.g., SG_REPR)
        
        Notes
        -----
        The clustering utility is populated by methods like `make_clustersYN` 
        and `make_graphYN`, which perform connected component analysis on 
        node subsets defined by eigenvector values.
        
        This structure is accessible via the `gcl` property and allows 
        efficient caching and retrieval of clustering results for different 
        eigenvectors and partitioning schemes.
        
        Examples
        --------
        The structure after clustering might look like::

            sg.gcl['eigV0']['+1'][SG_REPR]
            # Returns: (graphY, graphN)  - Tuple of subgraphs for nodes with eigV0 == +1
        """
        self.graph_clustering_utility = NestedDict()
    #
    def __init_paths__(
            self,
            path_data: Optional[Path] = None, 
            path_plot: Optional[Path] = None,
            make_dir_tree: bool = True,
            exist_ok: bool = True
    ) -> None:
        """
        Initialize paths for data, plots, and various subdirectories used by the
        SignedGraph class. If no paths are provided, default paths are used.

        Parameters
        ----------
        path_data : Path, optional
            The base path for data storage. Defaults to the global PATHDATA.
        path_plot : Path, optional
            The base path for plot storage. Defaults to the global PATHPLOT.
        make_dir_tree : bool, optional
            If True, creates the necessary directory structure for graph data.
        exist_ok : bool, optional
            If True, existing directories will not raise an error.

        Returns
        -------
        dict
            A dictionary containing the initialized paths.
        
        Notes
        -----
        This method creates directories only for graph-related data (graph, lrgsg,
        phtra, spect). Dynamics-specific directories (ising, voter, contact) are
        defined but not created until their respective dynamics classes 
        (IsingDynamics, VoterModel, ContactProcess) are instantiated.
        
        Ensures syshapePth exists. Many topology subclasses set 
        `self.syshapePth` during their own initialization. When SignedGraph 
        is created directly from a networkx.Graph (without those wrappers) 
        this attribute may be missing and would raise an AttributeError. 
        Provides sensible fallbacks: prefer an existing `self.syshape` 
        integer, otherwise use the number of nodes if `self.G` is already 
        present, or a generic placeholder.
        """
        #
        self.path_data = path_data or PATHDATA
        self.path_plot = path_plot or PATHPLOT
        self.path_sgdata = self.path_data / Path(self.sgpathn)
        #
        if not hasattr(self, 'syshapePth') or self.syshapePth is None:
            if hasattr(self, 'syshape') and isinstance(self.syshape, int):
                self.syshapePth = f"N={self.syshape}"
            else:
                try:
                    # If G was passed in, it might be available here
                    if self.G is not None:
                        self.syshapePth = f"N={len(self.G)}"
                    else:
                        self.syshapePth = "N=unknown"
                except Exception:
                    self.syshapePth = "N=unknown"

        # Create paths for graph-related subdirectories only
        # Dynamics-specific paths (ising, voter, contact) are created 
        # by their respective dynamics classes when instantiated
        self.subpath_list = []
        for p in PATHN_GRAPH_LIST:
            pfname = Path(p, self.syshapePth)
            setattr(self, f"path_{p}", self.path_sgdata / pfname)
            self.subpath_list.append(getattr(self, f"path_{p}"))
        
        # Initialize dynamics paths as None (will be set by dynamics classes)
        for p in PATHN_DYNAMICS_LIST:
            pfname = Path(p, self.syshapePth)
            setattr(self, f"path_{p}", self.path_sgdata / pfname)
        #
        if make_dir_tree: self.__init_dirs__(exist_ok)
    #
    def __init_reprdict__(self):
        """
        Initialize the graph representation dictionary.
        
        Populates the graph representation dictionary with available graph 
        representations and stores the list of representation keys. The 
        primary representation (on_g) is always included, while additional 
        representations from SG_LIST_REPR are added if they exist as 
        attributes.
        
        Notes
        -----
        Sets `self.gr` with the primary graph representation and any 
        additional representations found in SG_LIST_REPR. Also creates 
        `self.graph_reprs` as a list of all available representation keys.
        """
        self.gr[self.on_g] = getattr(self, self.on_g)
        for grd in SG_LIST_REPR:
            try:
                self.gr[grd] = getattr(self, grd)
            except AttributeError:
                pass
        self.graph_reprs = list(self.gr.keys())
    #
    def __init_weights__(
        self, 
        values: Union[float, dict] = 1., 
        on_g: Optional[str] = None
    ) -> None:
        """
        Initialize edge weights for the graph.
        
        Sets the 'weight' attribute for all edges in the specified graph
        representation and updates all graph representations to maintain
        consistency.
        
        Parameters
        ----------
        values : Union[float, dict], default 1.
            Weight value(s) to assign to edges. Can be either:
            - A single float value applied to all edges
            - A dictionary mapping edge tuples to weight values
        on_g : str, optional
            The graph representation to initialize weights on. If None,
            uses the primary graph representation (self.on_g).
            
        Notes
        -----
        After setting edge weights, this method calls `upd_GraphRepr_All`
        to ensure all graph representations remain synchronized.
        """
        on_g = on_g or self.on_g
        nx.set_edge_attributes(self.gr[on_g], values, 'weight')  # type: ignore[call-overload]
        self.upd_GraphRepr_All(on_g)
    #
    def __init_sgraph__(
            self, 
            init_weights_val: Union[float, dict] = 1.
    ) -> None:
        """
        Initialize the signed graph structure and edge sets.
        
        Sets up edge sets (all edges, flipped/negative edges, and positive 
        edges) based on whether the graph is loaded from file or newly 
        created. For loaded graphs, determines edge signs from existing 
        weights. For new graphs, randomly flips edges according to pflip 
        and optionally initializes weights.
        
        Parameters
        ----------
        init_weights_val : Union[float, dict], default 1.
            Initial weight value(s) for edges. Only used when creating a 
            new graph without pre-existing weights. Can be either:
            - A single float value applied to all edges
            - A dictionary mapping edge tuples to weight values
            
        Notes
        -----
        For loaded graphs (self.load_g=True):
        - Reads edge weights to identify negative edges (weight < 0)
        - Computes pflip from the ratio of negative to total edges
        - Populates fleset (negative edges) and lfeset (positive edges)
        
        For new graphs (self.load_g=False):
        - Randomly selects Ne_flips edges to flip according to pflip
        - Only initializes edge weights if none exist on the graph
        - Uses __init_weights__() to set initial weight values
        
        Always updates graph representations and matrices after initialization.
        """
        on_g = self.on_g
        self.eset[on_g] = set(list(self.gr[on_g].edges()))
        if self.load_g:
            edges_data = self.gr[on_g].edges(data=True)
            self.fleset[on_g] = set([
                (u, v) for u, v, _ in edges_data if _.get('weight', 1) < 0
            ])
            self.pflip = self.Ne_n / self.Ne
            self.Ne_flips = int(self.pflip * self.Ne)
            self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
        else:
            self.nflip = int(self.pflip * self.N)
            self.Ne_flips = int(self.pflip * self.Ne)
            self.fleset[on_g] = set(
                self.get_random_links(self.Ne_flips, on_g=on_g)
            )
            self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
            edges_with_weights = any(
                'weight' in data 
                for _, _, data in self.gr[on_g].edges(data=True)
            )
            if not edges_with_weights:
                self.__init_weights__(init_weights_val)
        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices()
    #
    def __init_loaded_graph__(
        self,
        path_data: Optional[Path] = None, 
        path_plot: Optional[Path] = None,
        on_g: str = SG_REPR
    ) -> None:
        """
        Initialize a SignedGraph from a previously saved/pickled graph object.
        
        This method is used to reinitialize a loaded graph object, setting up
        necessary attributes and paths after unpickling. It ensures backward
        compatibility with older pickled objects and properly initializes the
        graph representation dictionary and related structures.
        
        Parameters
        ----------
        path_data : Path, optional
            The base path for data storage. Defaults to the global PATHDATA.
        path_plot : Path, optional
            The base path for plot storage. Defaults to the global PATHPLOT.
        on_g : str, default SG_REPR
            The primary graph representation to use.
            
        Notes
        -----
        This method performs the following initialization steps:
        1. Resets the graph representation dictionary
        2. Sets load_g flag to True (indicating a loaded graph)
        3. Ensures backward compatibility with older pickled objects
        4. Initializes graph representations, signed graph structures, and paths
        
        This is typically called internally after unpickling a graph object.
        """
        self.graph_representation_dictionary = {}
        self.load_g = True
        self.on_g = on_g
        
        # Ensure backward compatibility with older pickled objects
        self.__ensure_required_attributes__()
        
        self.__init_reprdict__()
        self.__init_sgraph__()
        self.__init_paths__(path_data=path_data, path_plot=path_plot)
    #
    def _verify_pflip(self, pflip: float) -> None:
        """
        Verify if the provided pflip value is within the valid range.

        Parameters
        ----------
        pflip : float
            The probability of flipping edges in the signed graph.

        Raises
        ------
        ValueError
            If pflip is not within the range [LB_PFLIP, UB_PFLIP].
        """
        if not is_in_range(pflip, LB_PFLIP, UB_PFLIP):
            raise ValueError(SG_ERRMSG_PFLIP)
        else:
            self.pflip = pflip
            from ....config.funcs import peq_fstr
            self.peq_str = peq_fstr(pflip)
    #
    # load graph tools (imported from ._loaders)
    #
    from ._loaders import __load_graph__
    from ._loaders import _load_eigV
    from ._loaders import load_eigV_all
    #
    # export graph tools (imported from ._exports)
    #
    from ._exports import __export_graph__
    from ._exports import _export_edgel_bin
    from ._exports import export_adj_bin
    from ._exports import _export_eigV
    from ._exports import export_eigV_all
    from ._exports import export_ising_clust
    #
    # graph operations (imported from ._ongraph)
    #
    from ._ongraph import check_Ne_flips
    from ._ongraph import flip_sel_edges
    from ._ongraph import get_random_edges_from_set
    from ._ongraph import flip_random_fract_edges
    from ._ongraph import unflip_all
    from ._ongraph import set_edges_random_normal
    from ._ongraph import load_vec_on_nodes
    from ._ongraph import load_eigV_on_graph
    from ._ongraph import set_node_attributes
    #
    # graph representations (imported from ._representations)
    #
    from ._representations import upd_graph_matrices
    from ._representations import upd_edge_sets
    from ._representations import upd_Degree
    from ._representations import zip_reprNodes
    from ._representations import zip_reprEdges
    from ._representations import upd_NodeMap
    from ._representations import upd_EdgeMap
    from ._representations import upd_GraphRelabel
    from ._representations import upd_ReprMaps
    from ._representations import upd_GraphRepr_All
    #
    # topology tools (imported from ._topology)
    #
    from ._topology import get_laplacian
    from ._topology import get_signed_laplacian
    from ._topology import get_signed_laplacian_embedding
    from ._topology import make_rescaled_signed_laplacian
    from ._topology import nodes_in
    from ._topology import get_nodes_list
    from ._topology import get_node_attributes
    from ._topology import get_edge_data
    from ._topology import get_edge_mapping
    from ._topology import get_edge_color
    from ._topology import get_graph_neighbors
    from ._topology import get_adjacency_matrix
    from ._topology import get_degree_matrix
    from ._topology import get_abs_degree_matrix
    from ._topology import get_adjacency_matrix_for
    from ._topology import get_degree_matrix_for
    from ._topology import get_signed_degree_matrix_for
    from ._topology import get_laplacian_matrix_for
    from ._topology import get_signed_laplacian_matrix_for
    #
    # spectral tools (imported from ._spectral)
    #
    from ._spectral import get_eigV
    from ._spectral import get_eigV_check
    from ._spectral import get_eigV_binarized
    from ._spectral import get_eigV_bin_check
    from ._spectral import get_eigV_bin_check_list
    from ._spectral import get_sgspect_basis
    from ._spectral import compute_laplacian_spectrum
    from ._spectral import compute_laplacian_spectrum_weigV
    from ._spectral import compute_k_eigvV
    from ._spectral import compute_k_adj_eigvV
    from ._spectral import compute_adjacency_spectrum_weigV
    from ._spectral import make_eigV_transposed
    from ._spectral import make_eigV_column_major
    from ._spectral import make_adj_eigV_transposed
    from ._spectral import make_adj_eigV_column_major
    #
    # quantum propagator methods (imported from ._spectral)
    #
    from ._spectral import compute_quantum_propagator
    from ._spectral import quantum_walk_probabilities
    from ._spectral import quantum_observables_time_series
    #
    # information theory tools (imported from ._infotheory)
    #
    from ._infotheory import compute_signed_laplacian_entropy
    from ._infotheory import compute_renyi_entropy_profile
    from ._infotheory import get_entropy
    from ._infotheory import get_specific_heat
    from ._infotheory import get_entropy_derivative  # DEPRECATED
    from ._infotheory import get_renyi_results
    #
    # energy and dynamics tools (imported from ._dynamics)
    #
    from ._dynamics import compute_sksph_energy_eigV
    from ._dynamics import compute_sksph_energy_eigV_all
    from ._dynamics import get_sksph_energy_eigV
    from ._dynamics import get_all_sksph_energy_eigV
    from ._dynamics import compute_rbim_energy_eigV
    from ._dynamics import compute_rbim_energy_eigV_all
    from ._dynamics import get_rbim_energy_eigV
    from ._dynamics import get_all_rbim_energy_eigV
    #
    # partitioning tools (imported from ._partitioning)
    #
    from ._partitioning import get_subgraph_from_nodes
    from ._partitioning import get_nodes_subgraph_by_kv
    from ._partitioning import get_eigV_cluster_sizes
    from ._partitioning import get_cluster_distribution
    from ._partitioning import get_ferroAntiferro_regions
    from ._partitioning import make_graphYN
    from ._partitioning import make_clustersYN
    from ._partitioning import make_eigVclustersYN
    from ._partitioning import make_connected_component_by_edge
    from ._partitioning import handle_no_clust
    #
    # data cleaners and removers (imported from ._cleaners)
    #
    from ._cleaners import remove_ising_clust_files
    from ._cleaners import remove_edgl_file
    from ._cleaners import remove_adj_file
    from ._cleaners import remove_eigV_file
    from ._cleaners import remove_exported_files
    from ._cleaners import clean_gclutil
    #
    # order-parameter helpers (moved to _ordparams)
    from ._ordparams import compute_gap
    from ._ordparams import get_gap
    from ._ordparams import compute_pinf
    from ._ordparams import compute_gap_between

    #
    def get_p_fname(
        self, 
        who: str, 
        out_suffix: str = '', 
        ext: str = BIN, 
    ) -> str | Path:
        """Get the file name for exporting or importing graph data."""
        return build_p_fname(who, self.pflip, out_suffix=out_suffix, ext=ext)
    #
    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes without building the graph."""
        if hasattr(self, "syshape") and isinstance(self.syshape, int):
            return self.syshape
        raise NotImplementedError(
            "Subclasses must implement `get_expected_num_nodes`." )
    #
    def get_random_links(self, n: int = 1, only_in: str = '',
                         on_g: str = SG_REPR) -> list:
        """
        Get random edges from the graph.
        
        Parameters
        ----------
        n : int, default 1
            Number of random edges to return.
        only_in : str, default ''
            Filter for edge type: '', 'all' (all edges), '+', 'positive', '+1', 
            'plus' (positive edges only), or '-', 'negative', '-1', 'minus' 
            (negative edges only).
        on_g : str, default SG_REPR
            Graph representation to use.
            
        Returns
        -------
        list
            List of randomly selected edges.
        """
        match only_in:
            case ''|'all':
                return self.get_random_edges_from_set(n, self.eset[on_g], on_g)
            case '+'|'positive'|'+1'|'plus':
                return self.get_random_edges_from_set(n, self.lfeset[on_g], on_g)
            case '-'|'negative'|'-1'|'minus':
                return self.get_random_edges_from_set(n, self.fleset[on_g], on_g)
            case _:
                # Default to all edges if unknown filter
                return self.get_random_edges_from_set(n, self.eset[on_g], on_g)
    #
    # computations
    #
    #
    #
    #
    #
    # def compute_pinf(
    #         self,
    #         which: int = 0, 
    #         val: ConditionalPartitioningInput = 1,
    #         on_g: str = SG_REPR
    # ) -> None:
    #     """
    #     Compute the infinite cluster probability for a given eigenvector.
        
    #     Parameters
    #     ----------
    #     which : int, default 0
    #         Index of the eigenvector to analyze.
    #     val : Optional[ConditionalPartitioningInput], optional
    #         Conditional partitioning specification for clustering. Can be a 
    #         ConditionalPartitioning object, a number, a string (e.g., '>0'), 
    #         or a callable. If None, defaults to +1.
    #     on_g : str, default SG_REPR
    #         Graph representation to use.
    #     """
    #     # Default to +1 if no value provided
    #     if val is None:
    #         val = 1
        
    #     clustd = np.array(self.get_eigV_cluster_sizes(which, True, val, on_g))
    #     self.Pinf = clustd[0] / self.N
        
    #     # Avoid division by zero in variance calculation
    #     denominator = np.sum(clustd) - clustd[0]
    #     if denominator > 0:
    #         self.Pinf_var = np.sum(clustd@clustd - clustd[0]**2) / denominator
    #     else:
    #         self.Pinf_var = 0.0
        
    #     if hasattr(self, "Pinf_dict"):
    #         self.Pinf_dict[which] = (self.Pinf, self.Pinf_var)
    #     else:
    #         self.Pinf_dict = {which: (self.Pinf, self.Pinf_var)}
    # #
    # #
    # #
    # #
    
