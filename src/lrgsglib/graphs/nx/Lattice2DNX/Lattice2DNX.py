from os.path import join as pth_join
import warnings

import numpy as np
import networkx as nx
from networkx import convert_node_labels_to_integers, set_node_attributes, Graph

from ....config.const import *
from ....config.errwar import Lattice2DWarning
from ....utils.basic.iterables import compose
from ....utils.basic.numeric import is_positive_int
from ....utils.basic.arithmetic import adjust_to_even
from ..funcs import *
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX
from .generators_2d import *
from ..._shared._draw import draw as _draw_lattice2d
from ..._shared._nw_geometry import oriented_cell_edges, hub_central_edge
from ._nw_container import Lattice2DNXnwContainer
from ..._shared.animation.lattice2d import _Lattice2DAnimate, _Lattice2DPlot
#
class Lattice2DNX(SignedGraphNX):
    """
    Signed 2D lattice with optional periodic boundary conditions.

    The lattice is built in two representations:
    - ``H``: coordinate-labeled nodes (tuples)
    - ``G``: integer-labeled nodes (default representation)

    Parameters
    ----------
    side1 : int, default L2D_SIDE1
        Primary side length. If ``side2`` is provided, the larger of the two
        becomes ``side1``.
    geo : str, default L2D_GEO
        Lattice geometry. Accepts either a short alias or its full name
        (``L2D_GEO_SHRT_LIST`` / ``L2D_GEO_LIST``): ``"sqr"`` (squared, z=4),
        ``"tri"`` (triangular, z=6), ``"hex"`` (hexagonal, z=3),
        ``"oct_sqr"`` (octagonal_sqr / rhomb-octagonal, z=3),
        ``"kgm"`` (kagome, z=4), ``"tri_hex"`` (tri_hexagonal, z=3). The
        small-world variants ``"sqr_sw"`` / ``"tri_sw"`` are selected
        automatically when ``prew > 0``.
    side2 : int, default L2D_SIDE2
        Secondary side length. If omitted, a square lattice is used.
    pbc : bool, default L2D_PBC
        Whether to use periodic boundary conditions.
    fbc_val : float, default L2D_FBCV
        Fixed boundary value passed to lattice generators when applicable.
    stdFnameSFFX : str, default L2D_STDFN
        Filename suffix used in on-disk exports.
    sgpathn : str, default L2D_SGPATH
        Subpath for storing lattice data.
    with_positions : bool, default L2D_WITH_POS
        Whether to store node positions for plotting.
    bend_positions : bool, default L2D_BEND_POS
        Whether to bend positions for certain lattice geometries.
    prew : float, default L2D_PREW
        Rewiring probability for small-world variants.
    only_const_mode : bool, default L2D_ONLY_CONST_MODE
        If True, do not build the graph; only populate metadata.
    **kwargs : Any
        Forwarded to ``SignedGraph`` (e.g., ``pflip``, ``seed``,
        ``init_nw_dict``, ``path_data``, ``path_plot``).

    Notes
    -----
    ``pflip`` defines how many edges are marked for sign flips, but weights
    are not set negative until you call ``flip_random_fract_edges`` or
    ``flip_sel_edges``, unless weights already exist on the input graph.

    Examples
    --------
    >>> from lrgsglib.graphs.nx import Lattice2DNX
    >>> lat = Lattice2DNX(side1=10, geo="sqr", pflip=0.2, seed=1)
    >>> lat.flip_random_fract_edges()  # apply sign flips
    """
    eta_c = 1.128
    #
    def __init__(
        self,
        side1: int = L2D_SIDE1,
        geo: str = L2D_GEO,
        side2: int = L2D_SIDE2,
        pbc: bool = L2D_PBC,
        fbc_val: float = L2D_FBCV,
        stdFnameSFFX: str = L2D_STDFN,
        sgpathn: str = L2D_SGPATH,
        with_positions: bool = L2D_WITH_POS,
        bend_positions: bool = L2D_BEND_POS,
        prew: float = L2D_PREW,
        only_const_mode: bool = L2D_ONLY_CONST_MODE,
        **kwargs,
    ) -> None:
        self._verify_pflip(kwargs.get('pflip', 0.))
        self.only_const_mode = only_const_mode
        self.__init_side__(side1, side2)
        self.pbc = pbc
        self.prew = prew
        self.fbc_val = fbc_val
        self.__init_geo__(geo)
        #
        #
        self.__init_stdFname__(stdFnameSFFX)
        #
        self.sgpathn = pth_join(sgpathn, L2D_PATH_DICT[self.geo]) if sgpathn else L2D_PATH_DICT[self.geo]
        self.with_positions = with_positions
        self.bend_positions = bend_positions
        #
        if not only_const_mode:
            self.__init_lattice__()
        else:
            self.G = nx.Graph()
        super(Lattice2DNX, self).__init__(self.G, **kwargs)
    #
    def __init_side__(self, side1: int, side2: int) -> None:
        """
        Initialize the sides of the lattice, ensuring they are positive integers.

        Parameters
        ----------
        side1 : int
            The first side of the lattice.
        side2 : int
            The second side of the lattice.

        Raises
        ------
        ValueError
            If either side1 or side2 is not a positive integer.
        """
        if not is_positive_int(side1):
            raise ValueError("side1 must be a positive integer.")
        if side2 and not is_positive_int(side2):
            raise ValueError("side2 must be a positive integer.")
        if side2:
            self.side1, self.side2 = (side2, side1) if side2 > side1 else (side1, side2)
        else:
            self.side1 = side1
        #
    #
    def __init_geo__(self, geo: str) -> None:
        self.geo = geo
        if self.prew > 0.:
            self.geo = geo + '_sw'
        if geo not in L2D_GEO_LIST:
            if geo not in L2D_GEO_SHRT_LIST:
                warnings.warn(L2D_WARNMSG_GEO, Lattice2DWarning)
                self.geo = L2D_GEO
            else:
                self.geo = L2D_SHRT_GEO_DICT[self.geo]
        if not hasattr(self, 'side2'):
            if self.geo == 'hexagonal':
                self.side2 = self.side1
                self.side1 = adjust_to_even(self.side1/np.sqrt(3))
                if (self.side1 % 2 or self.side2 % 2) and self.pbc:
                    raise ValueError(L2D_ERRMSG_GEO)
            elif self.geo in ('kagome', 'tri_hexagonal'):
                self.side2 = self.side1
            else:
                self.side2 = self.side1
            # if (self.side1 % 2 or self.side2 % 2) and self.pbc:
            #     raise ValueError(DEFLattice2D_geoerrmsg)
                #
        # For rhomb_octagonal graphs each rhomb contributes four nodes
        self.node_multiplier = 4 if self.geo.startswith(L2D_SHRT_GEO_DICT['oct_sqr']) else 1
        total_nodes = self.node_multiplier * self.side1 * self.side2
        
        if self.side1 == self.side2:
            self.syshapePth = f"N={total_nodes}"
        elif self.side1 > self.side2:
            self.syshapePth = f"L1={self.side1}_L2={self.side2}"
        elif self.side2 > self.side1:
            self.syshapePth = f"L1={self.side2}_L2={self.side1}"
        if self.prew > 0.:
            self.syshapePth = self.syshapePth + f"_prew={self.prew:.3g}"
        
        # Set syshape to reflect actual number of nodes
        if self.node_multiplier > 1:
            # For oct_sqr: each side gets multiplied by 2 (since 2*2 = 4 nodes per rhomb)
            self.syshape = (2 * self.side1, 2 * self.side2)
        else:
            self.syshape = (self.side1, self.side2)
        #
        self.p_c = L2D_P_C_DICT.get(self.geo, float('nan'))
        if np.isfinite(self.p_c) and self.p_c > 0:
            self.r_c = np.sqrt(self.eta_c/(np.pi*self.p_c))
        else:
            self.r_c = float('nan')
    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.std_fname = L2D_GEO_SHRT_DICT[self.geo] + SFFX
    #
    def __init_lattice__(self) -> None:
        #
        if self.geo.startswith(L2D_SHRT_GEO_DICT['tri']):
            if self.prew == 0.:
                nxfunc = triangular_lattice_graph_FastPatch
            else:
                nxfunc = compose(triangular_lattice_graph_FastPatch, 
                                 rewire_edges_optimized, 
                                 g_kwargs={'prew': self.prew})
            self.z = 6
        elif self.geo.startswith(L2D_SHRT_GEO_DICT['sqr']):
            if self.prew == 0.:
                nxfunc = squared_lattice_graph_FastPatch
            else:
                nxfunc = compose(squared_lattice_graph_FastPatch, 
                                 rewire_edges_optimized, 
                                 g_kwargs={'prew': self.prew})
            self.z = 4
        elif self.geo == L2D_SHRT_GEO_DICT['hex']:
            self.z = 3
            nxfunc = hexagonal_lattice_graph_FastPatch
        elif self.geo.startswith(L2D_SHRT_GEO_DICT['oct_sqr']):
            self.z = 3
            nxfunc = rhomb_octagonal_graph_FastPatch
        elif self.geo == L2D_SHRT_GEO_DICT['kgm']:
            self.z = 4
            nxfunc = kagome_lattice_graph
        elif self.geo == L2D_SHRT_GEO_DICT['tri_hex']:
            self.z = 3
            nxfunc = tri_hexagonal_lattice_graph

        #
        self.H = nxfunc(
            self.side1, 
            self.side2, 
            periodic=self.pbc, 
            with_positions=self.with_positions,
            bend_positions=self.bend_positions
        )
        self.G = nx.convert_node_labels_to_integers(self.H)

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the lattice."""
        return int(self.node_multiplier * self.side1 * self.side2)
    #
    def degree_check(self, degree):
        return np.where(np.array(list(map(lambda x: x[1], 
                                          list(self.G.degree())))) != degree)
    #
    def _cell_posfn(self, on_g: str = L2D_ONREP):
        """Return a ``node -> (x, y)`` coordinate function for repr ``on_g``.

        ``sqr`` / ``tri`` / ``hex`` label nodes by integer ``(i, j)`` coordinates
        (directly for ``H``; via ``map_node`` for ``G``). ``oct_sqr`` / ``kgm`` /
        ``tri_hex`` use opaque labels (plain ints / edge-pair tuples) with no 2D
        coordinate, so the function returns ``None`` there and callers degrade to
        the geometry-free behaviour.
        """
        base = (lambda n: n) if on_g == "H" else self.map_node["H"]["G"].get

        def posfn(n):
            c = base(n)
            if (
                isinstance(c, tuple)
                and len(c) == 2
                and all(isinstance(x, (int, float)) for x in c)
            ):
                return (float(c[0]), float(c[1]))
            return None

        return posfn

    def get_central_edge(self, on_g: str = L2D_ONREP):
        """A central edge for the ``single*`` nwDict patterns.

        Uses the lattice's 2D coordinates when the geometry provides them
        (``sqr`` / ``tri`` / ``hex``): the edge of the centroid-nearest node whose
        midpoint is closest to the centroid — mirroring ``geometric_central_edge``
        on the GT side. Geometries whose labels carry no coordinate (``oct_sqr`` /
        ``kgm`` / ``tri_hex``) fall back to the geometry-free hub edge instead of
        raising (the old ``map_edge`` lookup ``KeyError``'d on them).
        """
        posfn = self._cell_posfn(on_g)
        pts = [(n, posfn(n)) for n in self.get_nodes_list(on_g)]
        pts = [(n, p) for (n, p) in pts if p is not None]
        if not pts:
            return hub_central_edge(self, on_g)
        arr = np.array([p for _, p in pts])
        centroid = arr.mean(axis=0)
        ci = int(np.argmin(((arr - centroid) ** 2).sum(axis=1)))
        central, cpos = pts[ci][0], arr[ci]
        nbrs = list(self.get_graph_neighbors(central, on_g))
        if not nbrs:
            return hub_central_edge(self, on_g)

        def _mid(j):
            pj = posfn(j)
            if pj is None:
                return np.inf
            return float(((0.5 * (np.array(pj) + cpos) - centroid) ** 2).sum())

        return (central, min(nbrs, key=_mid))

    def cell_edges(self, node, on_g: str = L2D_ONREP):
        """Coordinate-oriented ZERR cell on the **square** lattice (overrides the
        geometry-free base).

        On a square lattice the ``(i, j)`` coordinates are shared verbatim with
        ``Lattice2DGT``'s ``_pos`` grid, so picking the minimal face at the
        smallest CCW angle from ``+x`` (its ``+x/+y`` / north-east plaquette)
        gives distinct seeds distinct cells in the bulk *and* matches GT
        node-for-node. The other geometries either label nodes differently
        between engines (``tri``) or expose no shared 2D coordinate
        (``oct_sqr`` / ``kgm`` / ``tri_hex``, whose labels are opaque; ``hex``,
        whose sheared GT ``_pos`` disagrees with NX's index), so for them this
        defers to the engine-independent canonical lex-min face of the base.
        """
        if not str(self.geo).startswith("squared"):
            return super().cell_edges(node, on_g)
        box = None
        if self.pbc:
            box = (float(self.side1), float(self.side2))
        return oriented_cell_edges(self, node, on_g, self._cell_posfn(on_g), box)
    #
    nwContainer = Lattice2DNXnwContainer
    # Engine-agnostic 2D lattice drawing (shared with Lattice2DGT).
    # See lrgsglib.graphs._shared._draw.draw for the full signature.
    draw = _draw_lattice2d

    # Plotting / animation namespaces (shared with Lattice2DGT): the methods
    # live on accessors, so `lat.plot.draw(...)`, `lat.animate.states(...)`,
    # `lat.animate.largest_cluster(...)`, `lat.animate.make(...)`. See
    # lrgsglib.graphs._shared.animation.lattice2d for full signatures.
    plot = property(lambda self: _Lattice2DPlot(self))
    animate = property(lambda self: _Lattice2DAnimate(self))
