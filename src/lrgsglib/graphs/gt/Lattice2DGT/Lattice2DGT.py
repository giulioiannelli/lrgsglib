"""
Lattice2DGT: graph-tool implementation of 2D lattice signed graphs.

Mirrors the API of Lattice2DNX for easy switching between backends.
"""
from __future__ import annotations

from typing import Literal, Optional, Tuple, Union
import numpy as np

import graph_tool.all as gt

from ....config.const import (
    L2D_PBC,
    L2D_WITH_POS,
    L2D_ONREP,
    L2D_GEO_SHRT_LIST,
    L2D_SHRT_GEO_DICT,
    L2D_GEO_SHRT_DICT,
    L2D_ERRMSG_GEO,
    SG_INIT_NW_DICT,
)
from ....utils.basic.arithmetic import adjust_to_even
from ..SignedGraphGT import SignedGraphGT
from ..._shared._draw import draw as _draw_lattice2d
from . import _generators as _gen2d
from ._nw_container import Lattice2DGTnwContainer
from ..._shared._nw_container import geometric_central_edge


_SW_SUFFIX = "_sw"  # small-world variant marker (mirrors Lattice2DNX)

# All geometries are built through the same native-generator path so the GT and
# NX engines are interchangeable. Maps a base geo to its `_generators` builder;
# the lambda adapts to each generator's ``(side1, side2)`` / ``(n, m)`` order.
_GENERATORS = {
    'sqr': lambda s1, s2, per: _gen2d.squared(s1, s2, per),
    'tri': lambda s1, s2, per: _gen2d.triangular(s1, s2, per),
    'hex': lambda s1, s2, per: _gen2d.hexagonal(s1, s2, per),
    'oct_sqr': lambda s1, s2, per: _gen2d.rhomb_octagonal(s1, s2, per),
    'kgm': lambda s1, s2, per: _gen2d.kagome(s1, s2, per),
    'tri_hex': lambda s1, s2, per: _gen2d.tri_hexagonal(s1, s2, per),
}


class Lattice2DGT(SignedGraphGT):
    """
    2D lattice graph using graph-tool backend.

    Mirrors Lattice2DNX API for easy backend switching.

    Parameters
    ----------
    side1 : int
        Size in first dimension.
    side2 : int, optional
        Size in second dimension. If None, uses side1 (square lattice).
    geo : str
        Lattice geometry, as a short alias or its full name (mirrors
        Lattice2DNX): ``'sqr'`` (square, z=4), ``'tri'`` (triangular, z=6),
        ``'hex'`` (hexagonal, z=3), ``'oct_sqr'`` (rhomb-octagonal, z=3),
        ``'kgm'`` (kagome, z=4), ``'tri_hex'`` (tri-hexagonal, z=3).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    periodic : bool, default L2D_PBC
        Whether to use periodic boundary conditions. Defaults to the same
        value as Lattice2DNX (``L2D_PBC``) so both engines agree.
    prew : float
        Small-world rewiring probability (0.0 to 1.0). When ``prew > 0`` a
        fraction of edges is rewired, yielding the ``*_sw`` variant.
    seed : int, optional
        Random seed for reproducibility.

    Examples
    --------
    >>> lat = Lattice2DGT(side1=50, geo='sqr', pflip=0.3)
    >>> lat.flip_random_fract_edges()
    >>> print(f"Nodes: {lat.N}, Edges: {lat.num_edges}")

    >>> # Kagome lattice, identical topology to engine='nx'
    >>> lat_kgm = Lattice2DGT(side1=20, geo='kgm')

    Notes
    -----
    Every geometry is built by the NetworkX-free generators in
    ``_generators.py`` and matches the topology of the corresponding
    Lattice2DNX lattice exactly (verified by adjacency spectrum), so the two
    engines are interchangeable through the Lattice2D factory. The ``hex``
    geometry rescales the sides exactly as Lattice2DNX does (see
    ``_resolve_sides``).
    """

    # Supported base geometries (short aliases). The '_sw' small-world variants
    # are selected via ``prew > 0`` rather than passed directly.
    GEOMETRIES = frozenset(
        s for s in L2D_GEO_SHRT_LIST if not s.endswith(_SW_SUFFIX)
    )

    # Geometry-specific negative-link (nwDict) pattern container.
    nwContainer = Lattice2DGTnwContainer

    def __init__(
        self,
        side1: int,
        side2: Optional[int] = None,
        geo: str = 'sqr',
        pflip: float = 0.0,
        periodic: bool = L2D_PBC,
        prew: float = 0.0,
        seed: Optional[int] = None,
        with_positions: bool = L2D_WITH_POS,
        init_nw_dict: bool = SG_INIT_NW_DICT,
    ):
        # Normalise geo (accept short alias or full name) to a canonical base
        base, geo_label = self._normalise_geo(geo, prew)
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")
        if not 0.0 <= prew <= 1.0:
            raise ValueError(f"prew must be in [0, 1], got {prew}")

        self.geo = base
        self.periodic = periodic
        self.prew = prew
        self.with_positions = with_positions
        # Resolve side lengths exactly as Lattice2DNX does (swap + hex rescale)
        self.side1, self.side2 = self._resolve_sides(
            side1, side2, base, periodic)

        # Set seed (np.random drives the small-world rewiring, matching NX)
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Build the base graph through the shared native-generator path
        G, pos = self._assemble(
            *_GENERATORS[base](self.side1, self.side2, self.periodic))

        # Optional small-world rewiring (mirrors nx rewire_edges_optimized)
        if prew > 0.0:
            self._apply_small_world(G)

        # (Re)assign an all-positive sign property over every current edge
        G.edge_properties["sign"] = G.new_edge_property("int", val=1)

        # Store position property (always available for .draw()); also
        # register it as a graph vertex property so get_node_attributes('pos')
        # returns it, mirroring Lattice2DNX's with_positions semantics.
        self._pos = pos
        if with_positions:
            G.vertex_properties["pos"] = pos

        # Set syshapePth before super().__init__ so lazy path init uses it
        self._syshapePth = f"N={G.num_vertices()}"

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed,
                         sgpathn=f"l2d_{geo_label}_gt",
                         init_nw_dict=init_nw_dict)

    @classmethod
    def _normalise_geo(cls, geo: str, prew: float) -> Tuple[str, str]:
        """Map a short alias or full name to ``(base_geo, path_label)``.

        Accepts e.g. ``'sqr'`` / ``'squared'`` / ``'sqr_sw'`` / ``'squared_sw'``.
        An explicit ``_sw`` suffix or ``prew > 0`` flags the small-world variant
        used only for the on-disk path label.
        """
        if geo in L2D_SHRT_GEO_DICT:        # already a short alias
            short = geo
        elif geo in L2D_GEO_SHRT_DICT:      # full name -> short alias
            short = L2D_GEO_SHRT_DICT[geo]
        else:
            raise ValueError(
                f"geo must be one of {sorted(cls.GEOMETRIES)} (short alias) "
                f"or a full name in {L2D_GEO_SHRT_LIST}, got '{geo}'"
            )
        is_sw = short.endswith(_SW_SUFFIX)
        base = short[: -len(_SW_SUFFIX)] if is_sw else short
        if base not in cls.GEOMETRIES:
            raise ValueError(
                f"geo must be one of {sorted(cls.GEOMETRIES)}, got '{geo}'"
            )
        label = f"{base}{_SW_SUFFIX}" if (is_sw or prew > 0.0) else base
        return base, label

    def _resolve_sides(
        self, side1: int, side2: Optional[int], base: str, periodic: bool
    ) -> Tuple[int, int]:
        """Resolve ``(side1, side2)`` exactly as ``Lattice2DNX`` does.

        - if ``side2`` is given, the larger value becomes ``side1`` (a swap);
        - otherwise ``side2 = side1`` for every geometry except ``hex``, which
          rescales ``side1 -> adjust_to_even(side1 / sqrt(3))`` so the honeycomb
          is roughly isotropic. Periodic ``hex`` then requires even sides.
        """
        if side2 is not None:
            return (side2, side1) if side2 > side1 else (side1, side2)
        if base == 'hex':
            s2 = side1
            s1 = adjust_to_even(side1 / np.sqrt(3))
            if (s1 % 2 or s2 % 2) and periodic:
                raise ValueError(L2D_ERRMSG_GEO)
            return s1, s2
        return side1, side1

    def _assemble(
        self, n_nodes: int, edges, pos_array
    ) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Build a GT graph from an integer edge list and a position array.

        Used by the native exotic generators (oct_sqr, kgm, tri_hex). The sign
        property is (re)assigned by ``__init__`` after any rewiring.
        """
        G = gt.Graph(directed=False)
        G.add_vertex(n_nodes)
        G.add_edge_list([(int(u), int(v)) for u, v in edges])
        pos = G.new_vertex_property("vector<double>")
        for i in range(n_nodes):
            pos[G.vertex(i)] = [float(pos_array[i, 0]), float(pos_array[i, 1])]
        return G, pos

    def _apply_small_world(self, G: gt.Graph) -> None:
        """Rewire a fraction ``self.prew`` of edges in place (Watts-Strogatz).

        Mirrors ``lrgsglib.graphs.nx.funcs.rewire_edges_optimized``: one endpoint
        of each selected edge is reconnected to a random non-adjacent vertex,
        avoiding self-loops and duplicates. Driven by ``np.random`` so a set
        ``seed`` is reproducible.
        """
        n_v = G.num_vertices()
        base_edges = [(int(e.source()), int(e.target())) for e in G.edges()]
        existing = set()
        for u, v in base_edges:
            existing.add((u, v))
            existing.add((v, u))
        flags = np.random.rand(len(base_edges)) < self.prew
        for i, (u, v) in enumerate(base_edges):
            if not flags[i]:
                continue
            e = G.edge(u, v)
            if e is None:
                continue
            G.remove_edge(e)
            existing.discard((u, v))
            existing.discard((v, u))
            w = int(np.random.randint(0, n_v))
            attempts = 0
            while w == u or w == v or (u, w) in existing:
                w = int(np.random.randint(0, n_v))
                attempts += 1
                if attempts > n_v:
                    break
            if w != u and G.edge(u, w) is None:
                G.add_edge(u, w)
                existing.add((u, w))
                existing.add((w, u))

    @property
    def pos(self) -> gt.VertexPropertyMap:
        """Vertex position property map for visualization."""
        return self._pos

    def get_central_edge(self, on_g: str = L2D_ONREP) -> Tuple[int, int]:
        """Return a bulk edge nearest the geometric centre of the lattice.

        Mirrors the intent of ``Lattice2DNX.get_central_edge`` (a single central
        defect edge) in GT's integer node space. See
        ``graphs._shared._nw_container.geometric_central_edge``.
        """
        return geometric_central_edge(self, on_g)

    # Engine-agnostic 2D lattice drawing (shared with Lattice2DNX).
    # See lrgsglib.graphs._shared._draw.draw for the full signature.
    draw = _draw_lattice2d


__all__ = ["Lattice2DGT"]
