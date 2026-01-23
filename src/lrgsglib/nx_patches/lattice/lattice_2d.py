"""
Lattice2D: Backward-compatible 2D lattice wrapper.

This module provides the Lattice2D class that wraps the LatticeND base class
while maintaining full backward compatibility with the existing API.

Example
-------
>>> from lrgsglib.nx_patches.lattice import Lattice2D
>>> lat = Lattice2D(side1=10, geo='sqr', pflip=0.1, seed=42)
>>> lat.flip_random_fract_edges()
>>> lat.N
100
"""

from os.path import join as pth_join
from typing import Any
import random
import warnings

import networkx as nx
import numpy as np

from ...config.const import (
    L2D_BEND_POS,
    L2D_ERRMSG_GEO,
    L2D_FBCV,
    L2D_GEO,
    L2D_GEO_LIST,
    L2D_GEO_SHRT_DICT,
    L2D_GEO_SHRT_LIST,
    L2D_ONREP,
    L2D_ONLY_CONST_MODE,
    L2D_P_C_DICT,
    L2D_PATH_DICT,
    L2D_PBC,
    L2D_PREW,
    L2D_SGPATH,
    L2D_SHRT_GEO_DICT,
    L2D_SIDE1,
    L2D_SIDE2,
    L2D_STDFN,
    L2D_WARNMSG_GEO,
    L2D_WITH_POS,
    COUNT_XERR_PATTERNS,
)
from ...config.errwar import Lattice2DWarning
from ...utils.basic.arithmetic import adjust_to_even
from ...utils.basic.functions import compose
from ...utils.basic.numeric import is_positive_int
from ..funcs import get_neighbors_at_distance, rewire_edges_optimized

# Import generators
from ..Lattice2D.generators_2d import (
    hexagonal_lattice_graph_FastPatch,
    rhomb_octagonal_graph_FastPatch,
    squared_lattice_graph_FastPatch,
    triangular_lattice_graph_FastPatch,
)
from ..SignedGraph.SignedGraph import SignedGraph
from .lattice_nd import LatticeND, LatticeNwContainerBase


class Lattice2D(LatticeND):
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
        Lattice geometry. Common values: ``"sqr"``, ``"tri"``, ``"hex"``,
        ``"sqr_sw"``, ``"tri_sw"``, ``"oct_sqr"``.
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
    >>> from lrgsglib.nx_patches import Lattice2D
    >>> lat = Lattice2D(side1=10, geo="sqr", pflip=0.2, seed=1)
    >>> lat.flip_random_fract_edges()  # apply sign flips
    """

    eta_c = 1.128

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
        # Verify pflip before proceeding
        self._verify_pflip(kwargs.get('pflip', 0.0))

        # Store 2D-specific attributes before calling parent
        self.bend_positions = bend_positions
        self.prew = prew
        self.stdFnameSFFX = stdFnameSFFX

        # Process sides
        self._init_sides(side1, side2)

        # Process geometry (sets self.geo, node_multiplier, etc.)
        self._process_geometry(geo, pbc, sgpathn)

        # Initialize filename
        self._init_std_fname(stdFnameSFFX)

        # Call parent with computed dimensions
        super().__init__(
            dimensions=(self.side1, self.side2),
            geo=self.geo,
            pbc=pbc,
            fbc_val=fbc_val,
            with_positions=with_positions,
            only_const_mode=only_const_mode,
            **kwargs,
        )

    def _init_sides(self, side1: int, side2: int) -> None:
        """Initialize sides, ensuring positive integers and proper ordering."""
        if not is_positive_int(side1):
            raise ValueError("side1 must be a positive integer.")
        if side2 and not is_positive_int(side2):
            raise ValueError("side2 must be a positive integer.")

        if side2:
            self.side1, self.side2 = (
                (side2, side1) if side2 > side1 else (side1, side2)
            )
        else:
            self.side1 = side1
            self.side2 = side1

    def _process_geometry(self, geo: str, pbc: bool, sgpathn: str) -> None:
        """Process geometry string and set related attributes."""
        self.geo = geo

        # Handle small-world variants
        if self.prew > 0.0:
            self.geo = geo + '_sw'

        # Validate and normalize geometry
        if geo not in L2D_GEO_LIST:
            if geo not in L2D_GEO_SHRT_LIST:
                warnings.warn(L2D_WARNMSG_GEO, Lattice2DWarning)
                self.geo = L2D_GEO
            else:
                self.geo = L2D_SHRT_GEO_DICT[geo]

        # Handle hexagonal geometry special case
        if not hasattr(self, 'side2') or self.side2 == self.side1:
            if self.geo == 'hexagonal':
                self.side2 = self.side1
                self.side1 = adjust_to_even(self.side1 / np.sqrt(3))
                if (self.side1 % 2 or self.side2 % 2) and pbc:
                    raise ValueError(L2D_ERRMSG_GEO)
            elif not hasattr(self, 'side2'):
                self.side2 = self.side1

        # Set node multiplier for rhomb-octagonal
        self.node_multiplier = (
            4 if self.geo.startswith(L2D_SHRT_GEO_DICT['oct_sqr']) else 1
        )
        total_nodes = self.node_multiplier * self.side1 * self.side2

        # Set syshapePth
        if self.side1 == self.side2:
            self.syshapePth = f"N={total_nodes}"
        elif self.side1 > self.side2:
            self.syshapePth = f"L1={self.side1}_L2={self.side2}"
        else:
            self.syshapePth = f"L1={self.side2}_L2={self.side1}"

        if self.prew > 0.0:
            self.syshapePth = self.syshapePth + f"_prew={self.prew:.3g}"

        # Set syshape
        if self.node_multiplier > 1:
            self.syshape = (2 * self.side1, 2 * self.side2)
        else:
            self.syshape = (self.side1, self.side2)

        # Set physics constants
        self.p_c = L2D_P_C_DICT.get(self.geo, float('nan'))
        self.r_c = np.sqrt(self.eta_c / (np.pi * self.p_c))

        # Set path
        path_key = self.geo if self.geo in L2D_PATH_DICT else L2D_GEO
        self.sgpathn = (
            pth_join(sgpathn, L2D_PATH_DICT[path_key])
            if sgpathn
            else L2D_PATH_DICT[path_key]
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        geo_short = L2D_GEO_SHRT_DICT.get(self.geo, self.geo[:3])
        self.std_fname = geo_short + suffix

    def _get_geometry_registry(self) -> dict[str, dict]:
        """Return 2D geometry configuration registry."""
        return {
            'triangular': {'z': 6, 'gen': triangular_lattice_graph_FastPatch},
            'squared': {'z': 4, 'gen': squared_lattice_graph_FastPatch},
            'hexagonal': {'z': 3, 'gen': hexagonal_lattice_graph_FastPatch},
            'octagonal_sqr': {'z': 3, 'gen': rhomb_octagonal_graph_FastPatch},
        }

    def _init_geometry(self, geo: str) -> None:
        """Initialize geometry - already handled in _process_geometry."""
        pass  # Handled by _process_geometry for backward compatibility

    def _init_lattice(self) -> None:
        """Build the 2D lattice graph."""
        # Determine generator function
        if self.geo.startswith(L2D_SHRT_GEO_DICT['tri']):
            if self.prew == 0.0:
                nxfunc = triangular_lattice_graph_FastPatch
            else:
                nxfunc = compose(
                    triangular_lattice_graph_FastPatch,
                    rewire_edges_optimized,
                    g_kwargs={'prew': self.prew},
                )
            self.z = 6
        elif self.geo.startswith(L2D_SHRT_GEO_DICT['sqr']):
            if self.prew == 0.0:
                nxfunc = squared_lattice_graph_FastPatch
            else:
                nxfunc = compose(
                    squared_lattice_graph_FastPatch,
                    rewire_edges_optimized,
                    g_kwargs={'prew': self.prew},
                )
            self.z = 4
        elif self.geo == L2D_SHRT_GEO_DICT['hex']:
            self.z = 3
            nxfunc = hexagonal_lattice_graph_FastPatch
        elif self.geo.startswith(L2D_SHRT_GEO_DICT['oct_sqr']):
            self.z = 3
            nxfunc = rhomb_octagonal_graph_FastPatch
        else:
            # Default to squared
            nxfunc = squared_lattice_graph_FastPatch
            self.z = 4

        # Generate lattice
        self.H = nxfunc(
            self.side1,
            self.side2,
            periodic=self.pbc,
            with_positions=self.with_positions,
            bend_positions=self.bend_positions,
        )
        self.G = nx.convert_node_labels_to_integers(self.H)

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the lattice."""
        return int(self.node_multiplier * self.side1 * self.side2)

    def degree_check(self, degree: int) -> np.ndarray:
        """Return indices of nodes not matching the expected degree."""
        return np.where(
            np.array([d for _, d in self.G.degree()]) != degree
        )

    def get_central_edge(self, on_g: str = L2D_ONREP) -> tuple:
        """Return the central edge of the lattice."""
        cnode = (self.side1 // 2 - 1, self.side2 // 2)
        cnode_t = (self.side1 // 2, self.side2 // 2)

        if self.geo == 'triangular':
            cnode = (self.side2 // 2, self.side1 // 2 - 1)
            cnode_t = (self.side2 // 2, self.side1 // 2)

        edge_t = (cnode, cnode_t)

        if not self.H.has_edge(*edge_t):
            if self.geo == 'hexagonal':
                cnode = cnode_t
                cnode_t = (self.side1 // 2 + 1, self.side2 // 2)
                edge_t = (cnode, cnode_t)

        if on_g == 'G':
            return self.map_edge['G']['H'][edge_t]
        elif on_g == 'H':
            return edge_t

    class nwContainer(LatticeNwContainerBase):
        """Network container for 2D lattice patterns."""

        def __init__(
            self,
            l: SignedGraph,
            iterable: list = None,
            constant: Any = None,
            **kwargs,
        ):
            if iterable is None:
                iterable = []
            super().__init__(l, iterable, constant, **kwargs)

            # Add 2D-specific patterns
            self['singleZERR'] = {
                g: self.get_links_ZERR(self.centedge[g][0], g, self.l.geo)
                for g in self.rd
            }
            self['randZERR'] = {
                g: self.get_rand_pattern('ZERR', on_g=g)
                for g in self.rd
            }

        def get_links_ZERR(
            self, node: Any, on_g: str = L2D_ONREP, geometry: str = L2D_GEO
        ) -> list:
            """Get zero-error pattern based on geometry."""
            handlers = {
                'triangular': self.get_links_triangle,
                'squared': self.get_links_square,
                'hexagonal': self.get_links_hexagon,
            }
            return handlers.get(geometry, self.get_links_square)(node, on_g)

        def get_links_triangle(self, node: Any, on_g: str = L2D_ONREP) -> list:
            """Get triangular plaquette edges."""
            node2 = list(self.l.get_graph_neighbors(node, on_g))[0]
            common_neighbors = list(
                nx.common_neighbors(self.l.gr[on_g], node, node2)
            )
            try:
                node3 = common_neighbors[0]
                links = [(node, node2), (node2, node3), (node, node3)]
            except IndexError:
                links = [(node, node2)]
            return links

        def get_links_square(self, node: Any, on_g: str = L2D_ONREP) -> list:
            """Get square plaquette edges."""
            g = self.l.gr[on_g]
            neighbors = list(g.neighbors(node))

            for i in range(1, len(neighbors)):
                first_neighbor = neighbors[0]
                second_neighbor = neighbors[i]
                common_neighbors = (
                    set(g.neighbors(first_neighbor))
                    & set(g.neighbors(second_neighbor))
                )
                common_neighbors.discard(node)

                if common_neighbors:
                    common_neighbor = common_neighbors.pop()
                    return [
                        (node, first_neighbor),
                        (node, second_neighbor),
                        (first_neighbor, common_neighbor),
                        (second_neighbor, common_neighbor),
                    ]

            return [(node, neighbors[0]), (node, neighbors[1])] if len(neighbors) > 1 else [(node, neighbors[0])]

        def get_links_hexagon(self, node: int, on_g: str = L2D_ONREP) -> list:
            """Get hexagonal plaquette edges."""
            graph = self.l.gr[on_g]
            nodes_in_cycle = [node]
            node_nn = list(graph.neighbors(node))

            samp_node_nn_1 = node_nn[0]
            node_nn.remove(samp_node_nn_1)

            node_nn_1 = list(graph.neighbors(samp_node_nn_1))
            node_nn_1.remove(node)
            samp_node_nn_2 = node_nn_1[0]
            node_nn_1.remove(samp_node_nn_2)

            flag = True
            for nn in node_nn:
                if flag:
                    node_nn_1b = list(graph.neighbors(nn))
                    node_nn_1b.remove(node)
                    for i in node_nn_1b:
                        common_neighs = list(
                            nx.common_neighbors(graph, samp_node_nn_2, i)
                        )
                        if common_neighs:
                            nodes_in_cycle.extend([nn, i, common_neighs[0]])
                            flag = False
                            break

            try:
                nodes_in_cycle.extend([samp_node_nn_1, samp_node_nn_2])
            except IndexError:
                pass

            subH = graph.subgraph(nodes_in_cycle)
            return [tuple(sorted(edge)) for edge in subH.edges()]

        def get_rand_pattern(self, mode: str, on_g: str = L2D_ONREP) -> list:
            """Get random pattern of given mode."""
            if mode == "XERR":
                return self._get_rand_xerr_pattern(on_g)
            elif mode == "ZERR":
                return self._get_rand_zerr_pattern(on_g)
            return []

        def _get_rand_zerr_pattern(self, on_g: str) -> list:
            """Get random ZERR pattern based on geometry."""
            geo = self.l.geo
            if geo == 'squared':
                return [
                    k for i in self.rNodeFlip[on_g]
                    for k in self.get_links_square(i, on_g)
                ]
            elif geo == 'hexagonal':
                return [
                    k for i in self.rNodeFlip[on_g]
                    for k in self.get_links_hexagon(i, on_g)
                ]
            elif geo == 'triangular':
                return [
                    k for i in self.rNodeFlip[on_g]
                    for k in self.get_links_triangle(i, on_g)
                ]
            return []

        def get_links_rball(
            self,
            R: int = 1,
            center: Any = None,
            on_g: str = L2D_ONREP,
        ) -> set:
            """Get edges within radius R of center node."""
            graph = self.l.gr[on_g]
            if not center:
                center = self.centedge[on_g][0]
            neighs_to_flip = get_neighbors_at_distance(graph, center, R)
            return {
                (node, neighbor)
                for node in neighs_to_flip
                for neighbor in graph.neighbors(node)
            }

    def make_animation(
        self,
        fig,
        ax,
        frames,
        *,
        interval_ms: int = 50,
        cmap: str = "viridis",
        add_colorbar: bool = True,
        autoscale: bool = False,
        vmin: float | None = None,
        vmax: float | None = None,
        blit: bool = False,
    ):
        """Create animation from frames data."""
        from ..Lattice2D._animations import make_lattice2d_animation

        return make_lattice2d_animation(
            self,
            fig,
            ax,
            frames,
            interval_ms=interval_ms,
            cmap=cmap,
            add_colorbar=add_colorbar,
            autoscale=autoscale,
            vmin=vmin,
            vmax=vmax,
            blit=blit,
        )
