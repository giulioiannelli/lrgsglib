"""
Lattice3D: Backward-compatible 3D lattice wrapper.

This module provides the Lattice3D class that wraps the LatticeND base class
while maintaining full backward compatibility with the existing API.

Example
-------
>>> from lrgsglib.nx_patches.lattice import Lattice3D
>>> lat = Lattice3D(dim=8, geo='sc', pflip=0.1, seed=42)
>>> lat.flip_random_fract_edges()
>>> lat.N
512
"""

from os.path import join as pth_join
from typing import Any, Union, Tuple
import random
import warnings

import networkx as nx
from networkx import convert_node_labels_to_integers, set_node_attributes, Graph
import numpy as np

from ...config.const import (
    L3D_DIM,
    L3D_FBCV,
    L3D_GEO,
    L3D_GEO_BCC,
    L3D_GEO_DICT,
    L3D_GEO_FCC,
    L3D_GEO_LIST,
    L3D_GEO_SC,
    L3D_GEO_SHRT_LIST,
    L3D_ONREP,
    L3D_ONLY_CONST_MODE,
    L3D_PATH_DICT,
    L3D_PBC,
    L3D_PDIL,
    L3D_PHI,
    L3D_SGPATH,
    L3D_SHRT_GEO_DICT,
    L3D_STDFN,
    L3D_THETA,
    L3D_WITH_POS,
    L3D_WARNMSG_GEO,
    COUNT_XERR_PATTERNS,
)
from ...config.errwar import Lattice2DWarning
from ...utils.basic.functions import compose
from ...utils.basic.geometry import project_3d_to_2d
from ...utils.basic.numeric import is_positive_int
from ..funcs import LatticeND_graph_FastPatch, remove_edges

# Import 3D generators
from ..Lattice3D.generators_3d import generate_bcc_lattice, generate_fcc_lattice
from ..SignedGraph.SignedGraph import SignedGraph
from .lattice_nd import LatticeND, LatticeNwContainerBase


class Lattice3D(LatticeND):
    """
    Signed 3D lattice with optional periodic boundary conditions.

    The lattice is built in two representations:
    - ``H``: coordinate-labeled nodes (tuples)
    - ``G``: integer-labeled nodes (default representation)

    Parameters
    ----------
    dim : int | tuple[int, int, int], default L3D_DIM
        Lattice dimensions. If an int is provided, it is used for all axes.
    geo : str, default L3D_GEO
        Lattice geometry. Common values: ``"sc"``, ``"bcc"``, ``"fcc"``
        or their full names (``"simple_cubic"``, ``"body_centered"``,
        ``"face_centered"``).
    pbc : bool, default L3D_PBC
        Whether to use periodic boundary conditions.
    fbc_val : float, default L3D_FBCV
        Fixed boundary value passed to lattice generators when applicable.
    stdFnameSFFX : str, default L3D_STDFN
        Filename suffix used in on-disk exports.
    sgpathn : str, default L3D_SGPATH
        Subpath for storing lattice data.
    with_positions : bool, default L3D_WITH_POS
        Whether to store projected node positions for plotting.
    pdil : float, default L3D_PDIL
        Edge dilution probability for the simple cubic lattice.
    theta : float, default L3D_THETA
        Polar angle for 3D -> 2D projection when ``with_positions`` is True.
    phi : float, default L3D_PHI
        Azimuthal angle for 3D -> 2D projection when ``with_positions`` is True.
    only_const_mode : bool, default L3D_ONLY_CONST_MODE
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
    >>> from lrgsglib.nx_patches import Lattice3D
    >>> lat = Lattice3D(dim=6, geo="sc", pflip=0.1, seed=7)
    >>> lat.flip_random_fract_edges()  # apply sign flips
    """

    def __init__(
        self,
        dim: Union[int, Tuple[int, int, int]] = L3D_DIM,
        geo: str = L3D_GEO,
        pbc: bool = L3D_PBC,
        fbc_val: float = L3D_FBCV,
        stdFnameSFFX: str = L3D_STDFN,
        sgpathn: str = L3D_SGPATH,
        with_positions: bool = L3D_WITH_POS,
        pdil: float = L3D_PDIL,
        theta: float = L3D_THETA,
        phi: float = L3D_PHI,
        only_const_mode: bool = L3D_ONLY_CONST_MODE,
        **kwargs,
    ) -> None:
        # Store 3D-specific attributes
        self.pdil = pdil
        self.theta = theta
        self.phi = phi
        self.stdFnameSFFX = stdFnameSFFX

        # Process dimensions
        self._init_dim(dim)

        # Process geometry
        self._process_geometry(geo, sgpathn)

        # Initialize filename
        self._init_std_fname(stdFnameSFFX)

        # Call parent with computed dimensions
        super().__init__(
            dimensions=self.dim,
            geo=self.geo,
            pbc=pbc,
            fbc_val=fbc_val,
            with_positions=with_positions,
            only_const_mode=only_const_mode,
            **kwargs,
        )

        # Set positions after SignedGraph initialization
        if not only_const_mode and self.with_positions:
            self._set_positions()

    def _init_dim(self, dim: Union[int, Tuple[int, int, int]]) -> None:
        """Initialize dimensions from int or tuple."""
        if is_positive_int(dim):
            self.dim = (dim, dim, dim)
        elif (
            isinstance(dim, tuple)
            and len(dim) == 3
            and all(is_positive_int(x) for x in dim)
        ):
            self.dim = tuple(sorted(dim, reverse=True))
        else:
            raise ValueError(
                "dim must be a positive integer or a tuple of 3 positive integers"
            )

        self.dimL = list(self.dim)

    def _process_geometry(self, geo: str, sgpathn: str) -> None:
        """Process geometry string and set related attributes."""
        # Normalize geometry name
        self.geo = L3D_GEO_DICT.get(geo, geo)

        if self.geo not in L3D_GEO_LIST:
            if geo not in L3D_GEO_SHRT_LIST:
                warnings.warn(L3D_WARNMSG_GEO, Lattice2DWarning)
                self.geo = L3D_GEO_SC
            else:
                self.geo = L3D_SHRT_GEO_DICT.get(geo, L3D_GEO_SC)

        # Set node multiplier based on geometry
        if self.geo == L3D_GEO_BCC:
            self.node_multiplier = 2
        elif self.geo == L3D_GEO_FCC:
            self.node_multiplier = 4
        else:
            self.node_multiplier = 1

        # Set syshape and syshapePth
        self.syshape = self.dim
        total_nodes = self.node_multiplier * int(np.prod(self.dim))

        if all(x == self.dim[0] for x in self.dim):
            self.syshapePth = f"N={total_nodes}"
        else:
            dim_part = '_'.join([f"L{i}={side}" for i, side in enumerate(self.dim)])
            self.syshapePth = f"{dim_part}_N={total_nodes}"

        # Set path
        path_key = self.geo if self.geo in L3D_PATH_DICT else L3D_GEO_SC
        self.sgpathn = (
            pth_join(sgpathn, L3D_PATH_DICT[path_key])
            if sgpathn
            else L3D_PATH_DICT[path_key]
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        """Initialize standard filename."""
        self.std_fname = L3D_STDFN + suffix

    def _get_geometry_registry(self) -> dict[str, dict]:
        """Return 3D geometry configuration registry."""
        return {
            'simple_cubic': {'z': 6, 'gen': LatticeND_graph_FastPatch, 'node_mult': 1},
            'body_centered': {'z': 8, 'gen': generate_bcc_lattice, 'node_mult': 2},
            'face_centered': {'z': 12, 'gen': generate_fcc_lattice, 'node_mult': 4},
        }

    def _init_geometry(self, geo: str) -> None:
        """Initialize geometry - already handled in _process_geometry."""
        pass  # Handled by _process_geometry for backward compatibility

    def _init_lattice(self) -> None:
        """Build the 3D lattice graph."""
        if self.geo == L3D_GEO_SC:
            self.node_multiplier = 1
            if self.pdil == 0.0:
                nxfunc = LatticeND_graph_FastPatch
            else:
                nxfunc = compose(
                    LatticeND_graph_FastPatch,
                    remove_edges,
                    g_kwargs={'pdil': self.pdil},
                )
        elif self.geo == L3D_GEO_BCC:
            self.node_multiplier = 2
            nxfunc = generate_bcc_lattice
        elif self.geo == L3D_GEO_FCC:
            self.node_multiplier = 4
            nxfunc = generate_fcc_lattice
        else:
            raise ValueError(f"Unsupported geometry '{self.geo}'.")

        # Update metadata
        self.syshape = self.dim
        total_nodes = self.node_multiplier * int(np.prod(self.dim))

        if all(x == self.dim[0] for x in self.dim):
            self.syshapePth = f"N={total_nodes}"
        else:
            dim_part = '_'.join([f"L{i}={side}" for i, side in enumerate(self.dim)])
            self.syshapePth = f"{dim_part}_N={total_nodes}"

        # Generate lattice
        self.H = nxfunc(self.dim, periodic=self.pbc)
        self.G = convert_node_labels_to_integers(self.H)

    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes for the 3D lattice."""
        return int(self.node_multiplier * np.prod(self.dim))

    def _set_positions(self) -> None:
        """Set 2D projected positions on both H and G graphs."""
        # Set positions on H (tuple-labeled graph)
        pos_H = {
            node: project_3d_to_2d(*node, self.theta, self.phi)
            for node in self.H.nodes()
        }
        set_node_attributes(self.H, pos_H, 'pos')

        # Also set positions on G using the node mapping
        if (
            hasattr(self, 'map_node')
            and 'G' in self.map_node
            and 'H' in self.map_node['G']
        ):
            node_mapping = self.map_node['G']['H']  # Maps H nodes -> G nodes
            pos_G = {
                node_mapping[h_node]: position
                for h_node, position in pos_H.items()
            }
            set_node_attributes(self.G, pos_G, 'pos')
            # Update graph representation dictionary
            self.gr['G'] = self.G
            self.gr['H'] = self.H

    def get_central_edge(self, on_g: str = L3D_ONREP) -> tuple:
        """Return the central edge of the lattice."""
        cnode = (
            self.dimL[0] // 2 - 1,
            self.dimL[1] // 2,
            self.dimL[2] // 2,
        )
        cnode_t = (
            self.dimL[0] // 2,
            self.dimL[1] // 2,
            self.dimL[2] // 2,
        )
        edge_t = (cnode, cnode_t)

        if on_g == 'H':
            return edge_t
        elif on_g == 'G':
            return self.map_edge['G']['H'][edge_t]

    class nwContainer(LatticeNwContainerBase):
        """Network container for 3D lattice patterns."""

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
            # 3D lattice uses the base patterns only (no ZERR patterns)
