"""
SignedGraphGT: graph-tool backend for signed graphs.

This module provides a graph-tool implementation of SignedGraph that
mirrors the graphs.nx.SignedGraphNX API while leveraging graph-tool's
C++ performance.

Example
-------
>>> from lrgsglib.graphs.gt import SignedGraphGT
>>> sg = SignedGraphGT(n=100)  # Create from scratch
>>> # Or convert from NX
>>> from lrgsglib.graphs.gt._converters import nx_to_gt
>>> sg_gt = SignedGraphGT.from_networkx(nx_signed_graph)
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, diags
from scipy.sparse.csgraph import connected_components as sp_connected_components
from scipy.sparse.linalg import eigsh

try:
    import graph_tool as gt
    from graph_tool import Graph
    from graph_tool.spectral import adjacency, laplacian
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object  # Placeholder for type hints

import networkx as nx

from .._converters import nx_to_gt, gt_to_nx, get_laplacian_matrix_gt
from ..._shared._nw_container import NwContainer
from ..._shared._nw_geometry import hub_central_edge, elementary_cell_edges
from ..._shared._disorder import (
    Disorder,
    CompositeDisorder,
    as_disorder,
    plan_composite_ops,
)
from .._topology_ops import apply_rewiring, apply_dilution
from ....config.const import (
    BIN, PATHDATA, SG_REPR,
    PATHN_GRAPH_LIST, PATHN_DYNAMICS_LIST,
    SG_LAPL_SIGNED, SG_LAPL_RW, SG_LAPL_SYM,
    SG_LAPL_TYPES, SG_LAPL_DEFAULT_TYPE, SG_LAPL_RW_IMAG_TOL,
    SG_INIT_NW_DICT, SG_PREW, SG_PDIL, SG_DIL_EXTRACT_GIANT,
    SG_DISORDER,
)
from ....config.funcs import build_p_fname


def _real_cast_sorted_gt(w, V=None):
    """Sort a (possibly complex) spectrum ascending by real part and drop Im.

    The random-walk Laplacian ``L_rw`` is non-symmetric but isospectral to the
    symmetric PSD ``L_sym``, so its eigenvalues are real up to numerical noise.
    Warns if ``max|Im(lambda)|`` exceeds ``SG_LAPL_RW_IMAG_TOL``.
    """
    imag = float(np.max(np.abs(np.imag(w)))) if np.iscomplexobj(w) else 0.0
    if imag > SG_LAPL_RW_IMAG_TOL:
        import warnings

        warnings.warn(
            f"Non-symmetric (rw) Laplacian eigenvalues have |Im|={imag:.2e} > "
            f"{SG_LAPL_RW_IMAG_TOL:.1e}; casting to real.",
            RuntimeWarning,
            stacklevel=2,
        )
    order = np.argsort(np.real(w))
    w_real = np.ascontiguousarray(np.real(w)[order])
    if V is None:
        return w_real
    return w_real, np.ascontiguousarray(np.real(V[:, order]))


class SignedGraphGT:
    """
    Signed graph using graph-tool backend.

    This class provides a graph-tool based implementation of signed graphs,
    offering significant performance improvements for large graphs while
    maintaining API compatibility with SignedGraphNX.

    Parameters
    ----------
    G : Graph, optional
        Existing graph-tool Graph to wrap.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping.
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    G : Graph
        The underlying graph-tool Graph.
    N : int
        Number of nodes.
    pflip : float
        Fraction of edges marked for flipping.

    Notes
    -----
    Key differences from SignedGraphNX:
    - Uses graph-tool property maps instead of NetworkX edge attributes
    - Spectral operations are faster via graph-tool's C++ backend
    - Memory-efficient for large graphs

    Examples
    --------
    >>> from lrgsglib.graphs.gt import SignedGraphGT
    >>> # Create and flip edges
    >>> sg = SignedGraphGT.from_complete(n=50, pflip=0.2, seed=42)
    >>> sg.flip_random_fract_edges()
    >>> print(f"Negative edges: {sg.count_negative_edges()}")
    """

    # Negative-link pattern container class. Defaults to the engine-neutral
    # NwContainer so every GT graph gets the geometry-free patterns (rand,
    # randXERR) for free; geometric subclasses (e.g. Lattice2DGT) override with
    # a container that also builds single*/*ZERR. Mirrors SignedGraphNX.
    nwContainer: Optional[type] = NwContainer

    def __init__(
        self,
        G: Optional["Graph"] = None,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        path_data: Optional[Path] = None,
        sgpathn: str = "signed_graph_gt",
        init_nw_dict: bool = SG_INIT_NW_DICT,
        prew: float = SG_PREW,
        pdil: float = SG_PDIL,
        extract_giant_component: bool = SG_DIL_EXTRACT_GIANT,
        disorder: Optional[Union[str, "Disorder"]] = SG_DISORDER,
    ) -> None:
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # Resolve the disorder spec (support x coupling-law); a Disorder's own
        # pflip wins, otherwise the top-level pflip feeds the str/default path.
        self.disorder = as_disorder(disorder, pflip)
        self._pflip = self.disorder.pflip if self.disorder is not None else pflip
        self._seed = seed
        self._rng = np.random.default_rng(seed)

        if G is not None:
            self.G = G
        else:
            self.G = Graph(directed=False)

        # A graph handed in already carrying NEGATIVE signs (e.g. converted from
        # a signed NX graph via nx_to_gt, or loaded from disk) is respected
        # as-is: do NOT layer construction-time random disorder on top. Freshly
        # generated subclass graphs arrive all-positive (sign property all +1),
        # so they auto-flip normally. Mirrors the NX ``load_g`` guard.
        _incoming_signed = (
            G is not None
            and "sign" in self.G.edge_properties
            and bool((self.G.edge_properties["sign"].a < 0).any())
        )

        # ---- Topology disorder (engine-level, before signs) ----------------
        # Rewiring then dilution are applied to every GT graph here, so any
        # subclass that forwards prew/pdil gets them without bespoke code.
        # Done before the sign property is created so signs stay consistent.
        self._prew = prew
        self._pdil = pdil
        if prew > 0.0:
            self.G = apply_rewiring(self.G, prew)
        if pdil > 0.0:
            self.G = apply_dilution(self.G, pdil, extract_giant_component)
            # Dilution may prune+relabel vertices (giant component); keep a
            # geometric subclass's ``_pos`` consistent with the new graph.
            if "pos" in self.G.vertex_properties:
                self._pos = self.G.vertex_properties["pos"]

        # Ensure sign property exists
        if "sign" not in self.G.edge_properties:
            self.G.edge_properties["sign"] = self.G.new_edge_property("int", val=1)

        # Track edges marked for flipping
        self._flip_edges = set()

        # Set once a distributional coupling law writes continuous edge weights:
        # the signed adjacency then reads the 'weight' (double) property instead
        # of the ±1 'sign' (int) property. Discrete sign-flip graphs keep using
        # 'sign', so all existing behavior is byte-for-byte unchanged.
        self._continuous_couplings = False

        # Lazy neighbor cache for dynamics performance
        self._neighbor_cache: dict[int, list[tuple[int, float]]] | None = None

        # ---- Path infrastructure for DynamicsGraphProtocol ----
        self._path_data = path_data
        self._sgpathn = sgpathn
        # Paths are initialised lazily by _ensure_paths() on first access.
        self._paths_initialised = False

        # ---- Negative-link patterns (nwDict), built eagerly when requested ----
        # Mirrors SignedGraphNX: subclasses set the ``nwContainer`` class attr to
        # their geometry-specific container; here we instantiate it once.
        self.init_nw_dict = init_nw_dict
        need_nw = self.init_nw_dict or (
            self.disorder is not None and self.disorder.is_structured
        )
        if need_nw and self.nwContainer is not None:
            self.build_nw_dict(**self._nw_build_flags())

        # Realize the disorder at construction (auto-flip / distributional), so
        # SomeGraph(pflip=0.2) is signed on GT exactly as on NX. Uniform across
        # all GT subclasses now; the old per-subclass auto-flip calls are gone.
        # Skipped for graphs handed in pre-signed (respect the caller's signs).
        if self.disorder is not None and not _incoming_signed:
            self._apply_disorder(self.disorder)

    def get_central_edge(self, on_g: str = SG_REPR) -> Tuple[int, int]:
        """A central edge for the ``single*`` nwDict patterns.

        Base (geometry-free) implementation: the graph's hub — its
        highest-degree node and an edge to that node's highest-degree neighbour
        — so ``single`` / ``singleXERR`` / ``singleZERR`` supports work on any
        GT graph (ER/BA/Holme-Kim/…). Geometric subclasses (``Lattice2DGT`` /
        ``Lattice3DGT``) override with a coordinate-aware centre.
        """
        return hub_central_edge(self, on_g)

    def cell_edges(self, node, on_g: str = SG_REPR):
        """Edges of ``node``'s ZERR cell — the smallest cycle through it.

        Geometry-free baseline: the deterministic, engine-independent
        lexicographically-smallest minimal face (so graph-tool and NX agree).
        ``Lattice2DGT`` overrides with a coordinate-oriented plaquette so that
        distinct seed nodes own distinct cells in the bulk. Returns ``[]`` if
        ``node`` lies on no cycle.
        """
        return elementary_cell_edges(self, node, on_g)

    def _nw_build_flags(self) -> dict:
        """Container build-flag overrides for the requested disorder support.

        Geometry-free graphs default to a container that builds only
        ``rand``/``randXERR``; a structured ``single*`` / ``*ZERR`` support
        needs ``build_single`` / ``build_zerr`` forced on so the pattern exists.
        Only forwarded to engine-neutral :class:`NwContainer` subclasses.
        """
        if (
            self.disorder is None
            or not self.disorder.is_structured
            or not (
                isinstance(self.nwContainer, type)
                and issubclass(self.nwContainer, NwContainer)
            )
        ):
            return {}
        bs, bz = self.disorder.build_flags()
        flags: dict = {}
        if bs:
            flags["build_single"] = True
        if bz:
            flags["build_zerr"] = True
        return flags

    def build_nw_dict(self, **build_flags) -> None:
        """(Re)build ``self.nwDict`` from this graph's ``nwContainer``.

        Called automatically at construction when ``init_nw_dict=True``. Because
        the GT engine does **not** flip edges at construction (unlike NX —
        flips are applied later via ``flip_random_fract_edges``), the ``'rand'``
        pattern is a snapshot of the edges that are negative *at build time*.
        Call this again after flipping to refresh ``'rand'`` (the geometric
        ``single*`` / ``rand*`` star/cell patterns do not depend on signs).
        ``build_flags`` (``build_single`` / ``build_zerr``) force on-demand
        construction of the structured patterns for geometry-free graphs.
        """
        if self.nwContainer is None:
            raise NotImplementedError(
                f"{type(self).__name__} defines no nwContainer; "
                "nwDict patterns are unavailable for this graph type."
            )
        self.nwDict = self.nwContainer(self, **build_flags)

    def _disorder_support_edges(self, d: "Disorder") -> list:
        """Resolve a :class:`Disorder` support to a list of GT ``Edge`` objects.

        ``rand`` -> a random ``pflip`` fraction (sampled via
        :meth:`mark_flip_edges`); ``all`` -> every edge; a structured support
        (``randXERR``/…) -> the pattern stored in ``nwDict`` (built in
        ``__init__`` when the disorder is structured).
        """
        support = d.support
        if support == "rand":
            self.mark_flip_edges(d.pflip)
            return list(self._flip_edges)
        if support == "all":
            return list(self.G.edges())
        if d.is_registered_support:
            return self._edges_from_tuples(d.build_support(self, self._rng, SG_REPR))
        if d.is_structured:
            pattern = self.nwDict[support]
            tuples = pattern[SG_REPR] if isinstance(pattern, dict) else pattern
            return self._edges_from_tuples(tuples)
        return []

    def _edges_from_tuples(self, tuples) -> list:
        """Map ``(u, v)`` int tuples to existing GT ``Edge`` objects (skip absent)."""
        edges = []
        for (u, v) in tuples:
            e = self.G.edge(int(u), int(v))
            if e is not None:
                edges.append(e)
        return edges

    def _apply_disorder(self, d: "Disorder") -> None:
        """Realize a :class:`Disorder` on the GT graph at construction.

        ``flip`` SETs the support edges' ``sign`` to ``-1`` (idempotent), so
        ``SomeGraph(pflip=0.2)`` is signed on GT exactly as on NX. Distributional
        coupling laws require the ``sign`` property to be widened to ``double``
        (PLAN-Disorder-Model Stage 3) and are rejected here until then.
        """
        if d.is_none:
            return
        if isinstance(d, CompositeDisorder):
            self._apply_composite(d)
            return
        edges = self._disorder_support_edges(d)
        if not edges:
            return
        if d.is_flip:
            sign_prop = self.G.edge_properties["sign"]
            for e in edges:
                sign_prop[e] = -1
            # Record the realized set so a later no-arg flip_random_fract_edges()
            # re-asserts it idempotently (rather than re-sampling).
            self._flip_edges = set(edges)
        else:
            # Distributional couplings: write continuous values into a 'weight'
            # (double) property and switch matrix building onto it. The discrete
            # ±1 'sign' property is left intact; the signed Laplacian uses |J|
            # degrees, so continuous couplings need no further change.
            weight_prop = self._ensure_weight_prop()
            draws = d.draw(len(edges), self._rng)
            for e, w in zip(edges, draws):
                weight_prop[e] = float(w)
            self._continuous_couplings = True
        self.invalidate_cache()

    def _ensure_weight_prop(self):
        """The ``weight`` (double) edge property, created on first use (val 1.0)."""
        if "weight" not in self.G.edge_properties:
            self.G.edge_properties["weight"] = self.G.new_edge_property(
                "double", val=1.0
            )
        return self.G.edge_properties["weight"]

    def _apply_composite(self, comp: "CompositeDisorder") -> None:
        """Realize a :class:`CompositeDisorder` on GT.

        ``overlay`` applies each component in turn; ``mixture`` / the set-algebra
        modes reduce to combined write-ops via the engine-neutral planner. Because
        GT's signed matrix reads the ``weight`` property exclusively once any
        continuous law fires, a discrete ``flip`` op inside a composite that also
        carries a distributional component must ALSO write ``weight = -1``, else
        the flips would be read back as ``+1`` and silently lost.
        """
        if comp.mode == "overlay":
            for c in comp.components:
                self._apply_disorder(c)
            return
        force_weight = any(not c.is_flip for c in comp.components)
        ops = plan_composite_ops(comp, self, SG_REPR, self._resolve_edge_keys)
        for keys, law in ops:
            self._write_edge_keys(keys, law, force_weight=force_weight)

    def _resolve_edge_keys(self, d: "Disorder") -> set:
        """A component's support as canonical ``(min, max)`` int edge keys."""
        out = set()
        for e in self._disorder_support_edges(d):
            u, v = int(e.source()), int(e.target())
            out.add((u, v) if u <= v else (v, u))
        return out

    def _write_edge_keys(
        self, keys, law: "Disorder", force_weight: bool = False
    ) -> None:
        """Write one composite op (``flip`` or distributional) over ``keys``."""
        if not keys:
            return
        edges = self._edges_from_tuples(keys)
        if not edges:
            return
        if law.is_flip:
            sign_prop = self.G.edge_properties["sign"]
            for e in edges:
                sign_prop[e] = -1
            if force_weight:
                weight_prop = self._ensure_weight_prop()
                for e in edges:
                    weight_prop[e] = -1.0
            self._flip_edges = getattr(self, "_flip_edges", set()) | set(edges)
        else:
            weight_prop = self._ensure_weight_prop()
            draws = law.draw(len(edges), self._rng)
            for e, w in zip(edges, draws):
                weight_prop[e] = float(w)
            self._continuous_couplings = True
        self.invalidate_cache()

    def __getattr__(self, name: str) -> Any:
        """Trigger lazy path initialisation for ``path_*`` attributes."""
        if name.startswith("path_") and not name.startswith("_path_"):
            # Avoid recursion: check the flag directly on __dict__
            if not self.__dict__.get("_paths_initialised", True):
                self._ensure_paths()
                # After init, the attribute should exist now
                try:
                    return self.__dict__[name]
                except KeyError:
                    pass
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    @property
    def N(self) -> int:
        """Number of nodes."""
        return self.G.num_vertices()

    @property
    def num_edges(self) -> int:
        """Number of edges."""
        return self.G.num_edges()

    @property
    def Ne(self) -> int:
        """Number of edges (alias for num_edges, NX compatibility)."""
        return self.G.num_edges()

    @property
    def Ne_n(self) -> int:
        """Number of negative edges."""
        return self.count_negative_edges()

    @property
    def pflip(self) -> float:
        """Fraction of edges that are negative."""
        return self._pflip

    @pflip.setter
    def pflip(self, value: float) -> None:
        self._pflip = value

    @property
    def seed(self) -> Optional[int]:
        """Random seed used for graph generation."""
        return self._seed

    @seed.setter
    def seed(self, value: Optional[int]) -> None:
        self._seed = value

    @property
    def gr(self) -> Dict[str, "Graph"]:
        """Dictionary of graph representations (NX compatibility).

        Returns dict with 'G' key for primary graph. GT doesn't use
        multiple representations like NX.
        """
        return {"G": self.G}

    @property
    def graph_reprs(self) -> List[str]:
        """Available graph-representation keys (mirrors ``SignedGraphNX``).

        GT uses a single representation, so this is always ``["G"]``.
        """
        return list(self.gr.keys())

    @property
    def nflip(self) -> int:
        """Number of seed nodes used by the negative-link pattern containers.

        Mirrors ``SignedGraphNX.nflip`` (``int(pflip * N)``).
        """
        return int(self._pflip * self.N)

    @property
    def fleset(self) -> Dict[str, set]:
        """Flip-edge set per representation: currently-negative edges as
        ``(u, v)`` integer tuples (mirrors ``SignedGraphNX.fleset``).

        NX keeps ``fleset`` as a mutable dict updated on every flip; GT instead
        derives it from the live ``sign`` property so it always reflects the
        current state regardless of which flip path (``flip_random_fract_edges``,
        ``flip_sel_edges``, ``flip_edge``) produced the negatives.
        """
        sign_prop = self.G.edge_properties["sign"]
        neg = {
            (int(e.source()), int(e.target()))
            for e in self.G.edges()
            if sign_prop[e] == -1
        }
        return {"G": neg}

    @property
    def eigv(self) -> Optional[NDArray[np.floating]]:
        """Cached eigenvalues of signed Laplacian."""
        return getattr(self, "_eigv", None)

    @eigv.setter
    def eigv(self, value: NDArray[np.floating]) -> None:
        self._eigv = value

    @property
    def eigV(self) -> Optional[NDArray[np.floating]]:
        """Cached eigenvectors of signed Laplacian."""
        return getattr(self, "_eigV", None)

    @eigV.setter
    def eigV(self, value: NDArray[np.floating]) -> None:
        self._eigV = value

    @property
    def slp(self) -> Optional[NDArray[np.floating]]:
        """Signed Laplacian matrix (cached on first access)."""
        if not hasattr(self, "_slp") or self._slp is None:
            self._slp = self.get_signed_laplacian()
        return self._slp

    @property
    def adj(self) -> Optional[NDArray[np.floating]]:
        """Signed adjacency matrix (cached on first access)."""
        if not hasattr(self, "_adj") or self._adj is None:
            self._adj = self.get_signed_adjacency()
        return self._adj

    @classmethod
    def from_networkx(
        cls,
        G_nx: nx.Graph,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> "SignedGraphGT":
        """
        Create SignedGraphGT from NetworkX graph.

        Parameters
        ----------
        G_nx : nx.Graph
            NetworkX graph (with optional 'sign' edge attribute).
        pflip : float, default 0.0
            Fraction of edges to mark for flipping.
        seed : int, optional
            Random seed.

        Returns
        -------
        SignedGraphGT
            New instance wrapping converted graph.
        """
        G_gt = nx_to_gt(G_nx)
        return cls(G=G_gt, pflip=pflip, seed=seed)

    @classmethod
    def from_complete(
        cls,
        n: int,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> "SignedGraphGT":
        """
        Create complete (fully connected) signed graph.

        Parameters
        ----------
        n : int
            Number of nodes.
        pflip : float, default 0.0
            Fraction of edges to mark for flipping.
        seed : int, optional
            Random seed.

        Returns
        -------
        SignedGraphGT
            Complete graph with all positive edges.
        """
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")

        from graph_tool.generation import complete_graph
        G = complete_graph(n)
        return cls(G=G, pflip=pflip, seed=seed)

    @classmethod
    def from_lattice(
        cls,
        shape: tuple,
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> "SignedGraphGT":
        """
        Create lattice signed graph.

        Parameters
        ----------
        shape : tuple
            Shape of lattice (e.g., (10, 10) for 2D).
        periodic : bool, default True
            Whether to use periodic boundary conditions.
        pflip : float, default 0.0
            Fraction of edges to mark for flipping.
        seed : int, optional
            Random seed.

        Returns
        -------
        SignedGraphGT
            Lattice graph with all positive edges.
        """
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")

        from graph_tool.generation import lattice
        G = lattice(shape, periodic=periodic)
        return cls(G=G, pflip=pflip, seed=seed)

    def to_networkx(self) -> nx.Graph:
        """
        Convert to NetworkX graph.

        Returns
        -------
        nx.Graph
            NetworkX graph with 'sign' edge attribute.
        """
        return gt_to_nx(self.G)

    def mark_flip_edges(self, fraction: Optional[float] = None) -> None:
        """
        Mark a fraction of edges for sign flipping.

        Parameters
        ----------
        fraction : float, optional
            Fraction to mark. If None, uses self.pflip.
        """
        frac = fraction if fraction is not None else self.pflip
        if frac <= 0:
            return

        edges = list(self.G.edges())
        n_flip = int(len(edges) * frac)

        if n_flip > 0:
            indices = self._rng.choice(len(edges), size=n_flip, replace=False)
            self._flip_edges = {edges[i] for i in indices}

    def flip_random_fract_edges(self, fraction: Optional[float] = None) -> None:
        """
        Flip signs of a random fraction of edges.

        Parameters
        ----------
        fraction : float, optional
            Fraction to flip. If None, uses self.pflip.
        """
        # Idempotent for the no-arg case: re-assert the current flip set rather
        # than re-sampling, so a redundant call after construction-time
        # auto-flip is a no-op. An explicit ``fraction`` re-samples; an empty
        # flip set (deferred disorder=None workflow) samples on first call.
        if fraction is not None or not self._flip_edges:
            self.mark_flip_edges(fraction)
        sign_prop = self.G.edge_properties["sign"]

        for e in self._flip_edges:
            sign_prop[e] = -1

        # Invalidate cached matrices
        self.invalidate_cache()

    def flip_edge(self, u: int, v: int) -> None:
        """
        Flip sign of a specific edge.

        Parameters
        ----------
        u, v : int
            Vertex indices.
        """
        e = self.G.edge(u, v)
        if e is not None:
            self.G.edge_properties["sign"][e] *= -1
            self.invalidate_cache()

    def flip_sel_edges(
        self,
        edges: List[Tuple[int, int]],
    ) -> None:
        """
        Flip signs of selected edges.

        Parameters
        ----------
        edges : list of (u, v) tuples
            Edges to flip.
        """
        sign_prop = self.G.edge_properties["sign"]

        for u, v in edges:
            e = self.G.edge(u, v)
            if e is not None:
                sign_prop[e] = -sign_prop[e]

        self.invalidate_cache()

    def unflip_all(self) -> None:
        """Reset all edge signs to positive (+1)."""
        sign_prop = self.G.edge_properties["sign"]

        for e in self.G.edges():
            sign_prop[e] = 1

        self._flip_edges = set()
        self.invalidate_cache()

    def get_random_links(
        self,
        n: int = 1,
        only_in: Literal["", "all", "+", "positive", "-", "negative"] = "",
    ) -> List[Tuple[int, int]]:
        """
        Get random edges from the graph.

        Parameters
        ----------
        n : int, default 1
            Number of random edges to return.
        only_in : str, default ''
            Filter: '', 'all' (any edge), '+', 'positive' (positive only),
            '-', 'negative' (negative only).

        Returns
        -------
        list of (u, v) tuples
            Randomly selected edges.
        """
        sign_prop = self.G.edge_properties["sign"]

        if only_in in ("", "all"):
            edges = list(self.G.edges())
        elif only_in in ("+", "positive", "+1", "plus"):
            edges = [e for e in self.G.edges() if sign_prop[e] == 1]
        elif only_in in ("-", "negative", "-1", "minus"):
            edges = [e for e in self.G.edges() if sign_prop[e] == -1]
        else:
            edges = list(self.G.edges())

        if not edges:
            return []

        n = min(n, len(edges))
        indices = self._rng.choice(len(edges), size=n, replace=False)

        return [(int(edges[i].source()), int(edges[i].target())) for i in indices]

    def get_edge_sign(self, u: int, v: int) -> int:
        """
        Get sign of an edge.

        Parameters
        ----------
        u, v : int
            Vertex indices.

        Returns
        -------
        int
            Edge sign (+1 or -1), or 0 if edge doesn't exist.
        """
        e = self.G.edge(u, v)
        if e is not None:
            return self.G.edge_properties["sign"][e]
        return 0

    def count_negative_edges(self) -> int:
        """
        Count number of negative edges.

        Returns
        -------
        int
            Number of edges with sign = -1.
        """
        prop = self._coupling_edge_prop()
        return sum(1 for e in self.G.edges() if prop[e] < 0)

    def count_positive_edges(self) -> int:
        """
        Count number of positive edges.

        Returns
        -------
        int
            Number of edges with sign = +1.
        """
        return self.num_edges - self.count_negative_edges()

    def _coupling_edge_prop(self):
        """Edge property carrying the signed couplings for matrix building.

        Discrete sign-flip graphs use the ±1 ``'sign'`` (int) property; graphs
        with a distributional coupling law use the continuous ``'weight'``
        (double) property written by :meth:`_apply_disorder`. Gated on
        ``_continuous_couplings`` so existing behavior is unchanged.
        """
        if getattr(self, "_continuous_couplings", False) and (
            "weight" in self.G.edge_properties
        ):
            return self.G.edge_properties["weight"]
        return self.G.edge_properties["sign"]

    def get_signed_adjacency(self) -> np.ndarray:
        """
        Get signed adjacency matrix.

        Returns
        -------
        np.ndarray
            Signed adjacency matrix ``A`` where ``A[i,j]`` is the signed
            coupling of edge ``(i,j)`` (``±1`` by default, or the continuous
            coupling when a distributional disorder law is in effect).
        """
        return adjacency(self.G, weight=self._coupling_edge_prop()).toarray()

    def get_adjacency_matrix(self) -> np.ndarray:
        """
        Get unsigned adjacency matrix.

        Returns
        -------
        np.ndarray
            Adjacency matrix (1 for edges, 0 otherwise).
        """
        return adjacency(self.G).toarray()

    def get_signed_laplacian(self) -> np.ndarray:
        """
        Get signed Laplacian matrix.

        Returns
        -------
        np.ndarray
            Signed Laplacian L_s = D_s - A_signed where D_s[i,i] = sum_j |A_ij|.

        Notes
        -----
        The signed Laplacian uses absolute degree values to ensure
        proper spectral properties. This matches SignedGraphNX.
        """
        # Get signed adjacency
        A_signed = self.get_signed_adjacency()

        # Compute absolute degree (sum of absolute edge weights)
        degrees = np.sum(np.abs(A_signed), axis=1)

        # Build signed Laplacian: L_s = D_s - A_signed
        D_s = np.diag(degrees)
        return D_s - A_signed

    def get_laplacian(self) -> np.ndarray:
        """
        Get unsigned Laplacian matrix.

        Returns
        -------
        np.ndarray
            Standard Laplacian L = D - A.
        """
        return laplacian(self.G).toarray()

    def get_degree_matrix(self) -> np.ndarray:
        """
        Get degree matrix.

        Returns
        -------
        np.ndarray
            Diagonal matrix with node degrees.
        """
        degrees = np.array([v.out_degree() for v in self.G.vertices()])
        return np.diag(degrees)

    def get_abs_degree_matrix(self) -> np.ndarray:
        """
        Get absolute degree matrix (for signed Laplacian).

        Returns
        -------
        np.ndarray
            Diagonal matrix with sum of absolute edge weights.
        """
        A_signed = self.get_signed_adjacency()
        degrees = np.sum(np.abs(A_signed), axis=1)
        return np.diag(degrees)

    def get_signed_rw_laplacian(
        self, sym: bool = True, format: str = "dense"
    ) -> np.ndarray:
        """
        Random-walk normalized signed Laplacian (Kunegis 2010, §3.3/§3.4).

        Parameters
        ----------
        sym : bool, default True
            ``True``  -> ``L_sym = I - D_s^{-1/2} A D_s^{-1/2}`` (symmetric, PSD).
            ``False`` -> ``L_rw  = I - D_s^{-1} A`` (random walk, NOT symmetric).
        format : str, default 'dense'
            Accepted for NX-API symmetry but ignored (GT is always dense).

        Notes
        -----
        ``D_s[i,i] = sum_j |A_ij|`` is the signed degree. The two variants are
        isospectral, with real eigenvalues in ``[0, 2]``. Isolated nodes get a
        zero inverse (that row stays equal to the identity row, eigenvalue 1).
        """
        A = self.get_signed_adjacency()
        deg = np.abs(A).sum(axis=1)
        n = A.shape[0]
        eye = np.eye(n)
        if sym:
            inv = np.where(deg > 0, 1.0 / np.sqrt(deg), 0.0)
            return eye - (inv[:, None] * A) * inv[None, :]
        inv = np.where(deg > 0, 1.0 / deg, 0.0)
        return eye - inv[:, None] * A

    def _laplacian_operator(self, laplacian_type: str = SG_LAPL_DEFAULT_TYPE):
        """Resolve ``laplacian_type`` to ``(matrix, is_symmetric)``.

        ``'signed'`` -> combinatorial ``D_s - A`` (symmetric, default).
        ``'sym'``    -> symmetric normalized ``L_sym`` (symmetric).
        ``'rw'``     -> random-walk ``L_rw`` (NOT symmetric; general solver).
        """
        if laplacian_type == SG_LAPL_SIGNED:
            return self.get_signed_laplacian(), True
        if laplacian_type == SG_LAPL_SYM:
            return self.get_signed_rw_laplacian(sym=True), True
        if laplacian_type == SG_LAPL_RW:
            return self.get_signed_rw_laplacian(sym=False), False
        raise ValueError(
            f"laplacian_type must be one of {SG_LAPL_TYPES}, got {laplacian_type!r}"
        )

    def get_nodes_list(self, on_g: Optional[str] = None) -> List[int]:
        """
        Get list of node indices.

        Parameters
        ----------
        on_g : str, optional
            Accepted for cross-engine API parity (GT has a single
            representation); ignored.

        Returns
        -------
        list[int]
            Node indices [0, 1, ..., N-1].
        """
        return list(range(self.N))

    def get_edges_list(self) -> List[Tuple[int, int]]:
        """
        Get list of edges as (source, target) tuples.

        Returns
        -------
        list of tuples
            All edges in the graph.
        """
        return [(int(e.source()), int(e.target())) for e in self.G.edges()]

    def get_edge_signs(self) -> Dict[Tuple[int, int], int]:
        """
        Get dictionary mapping edges to their signs.

        Returns
        -------
        dict
            {(u, v): sign} for all edges.
        """
        sign_prop = self.G.edge_properties["sign"]
        return {
            (int(e.source()), int(e.target())): sign_prop[e]
            for e in self.G.edges()
        }

    def get_laplacian_spectrum(
        self,
        k: Optional[int] = None,
        backend: Literal["numpy", "scipy", "cupy"] = "numpy",
        laplacian_type: str = SG_LAPL_DEFAULT_TYPE,
    ) -> np.ndarray:
        """
        Compute Laplacian eigenvalues.

        Parameters
        ----------
        k : int, optional
            Number of eigenvalues to compute. If None, compute all.
        backend : str, default 'numpy'
            Computation backend ('numpy', 'scipy', 'cupy').
        laplacian_type : str, default 'signed'
            Which Laplacian to diagonalize ('signed', 'sym', 'rw').

        Returns
        -------
        np.ndarray
            Sorted eigenvalues.
        """
        L, is_sym = self._laplacian_operator(laplacian_type)

        if k is not None and k < self.N:
            L_sparse = csr_matrix(L)
            if is_sym:
                eigenvalues, _ = eigsh(L_sparse, k=k, which="SM")
            else:
                from scipy.sparse.linalg import eigs
                eigenvalues, _ = eigs(L_sparse, k=k, which="SR")
                eigenvalues = _real_cast_sorted_gt(eigenvalues)
        elif not is_sym:
            # Full non-symmetric (rw) spectrum -> dense general solver
            from scipy.linalg import eigvals as _sp_eigvals
            eigenvalues = _real_cast_sorted_gt(_sp_eigvals(L))
        elif backend == "cupy":
            try:
                import cupy as cp
                eigenvalues = cp.asnumpy(cp.linalg.eigvalsh(cp.asarray(L)))
            except ImportError:
                eigenvalues = np.linalg.eigvalsh(L)
        elif backend == "scipy":
            from scipy.linalg import eigvalsh
            eigenvalues = eigvalsh(L)
        else:
            eigenvalues = np.linalg.eigvalsh(L)

        return np.sort(eigenvalues)

    def compute_laplacian_spectrum_weigV(
        self,
        backend: Literal["numpy", "scipy", "cupy"] = "numpy",
        laplacian_type: str = SG_LAPL_DEFAULT_TYPE,
    ) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
        """
        Compute full eigendecomposition of the signed Laplacian.

        Parameters
        ----------
        backend : str, default 'numpy'
            Computation backend. 'numpy' uses np.linalg.eigh.
            'scipy' uses scipy.linalg.eigh.
            'cupy' uses GPU acceleration if available.
        laplacian_type : str, default 'signed'
            Which Laplacian to diagonalize ('signed', 'sym', 'rw'). For 'rw'
            (non-symmetric) a general solver is used and eigenvectors are the
            true (non-orthonormal) random-walk eigenvectors.

        Returns
        -------
        tuple[ndarray, ndarray]
            (eigenvalues, eigenvectors) where eigenvalues are sorted
            ascending and eigenvectors are column-major (eigV[:, i] is
            the i-th eigenvector).

        Notes
        -----
        Results are cached in self.eigv and self.eigV.
        """
        L, is_sym = self._laplacian_operator(laplacian_type)

        if not is_sym:
            # Non-symmetric (rw): general dense solver, real-cast & sorted.
            # cupy lacks a general (non-Hermitian) eig, so fall back to numpy.
            if backend == "scipy":
                from scipy.linalg import eig as _sp_eig
                eigenvalues, eigenvectors = _sp_eig(L)
            else:
                eigenvalues, eigenvectors = np.linalg.eig(L)
            self.eigv, self.eigV = _real_cast_sorted_gt(eigenvalues, eigenvectors)
            self._spectrum_laplacian_type = laplacian_type
            return self.eigv, self.eigV

        if backend == "numpy":
            eigenvalues, eigenvectors = np.linalg.eigh(L)
        elif backend == "scipy":
            from scipy.linalg import eigh

            eigenvalues, eigenvectors = eigh(L)
        elif backend == "cupy":
            try:
                import cupy as cp

                L_gpu = cp.asarray(L)
                eigenvalues, eigenvectors = cp.linalg.eigh(L_gpu)
                eigenvalues = cp.asnumpy(eigenvalues)
                eigenvectors = cp.asnumpy(eigenvectors)
            except ImportError:
                # Fallback to numpy
                eigenvalues, eigenvectors = np.linalg.eigh(L)
        else:
            raise ValueError(f"Unknown backend: {backend}")

        # Sort by eigenvalue (ascending)
        idx = np.argsort(eigenvalues)
        self.eigv = eigenvalues[idx]
        self.eigV = eigenvectors[:, idx]
        self._spectrum_laplacian_type = laplacian_type

        return self.eigv, self.eigV

    def compute_k_eigvV(
        self,
        k: int = 1,
        which: Literal["SM", "LM"] = "SM",
        laplacian_type: str = SG_LAPL_DEFAULT_TYPE,
        **kwargs: Any,
    ) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
        """
        Compute k smallest or largest eigenvalues/eigenvectors.

        Results are stored on ``self.eigv`` and ``self.eigV`` (column-major).

        Parameters
        ----------
        k : int, default 1
            Number of eigenvalues to compute.
        which : str, default 'SM'
            'SM' for smallest magnitude, 'LM' for largest magnitude.
        laplacian_type : str, default 'signed'
            Which Laplacian to diagonalize ('signed', 'sym', 'rw').
        **kwargs
            Ignored (accepted for NX API compatibility, e.g. ``typf``).

        Returns
        -------
        tuple[ndarray, ndarray]
            (eigenvalues, eigenvectors) with k values.
        """
        L, is_sym = self._laplacian_operator(laplacian_type)
        L_sparse = csr_matrix(L)

        if is_sym:
            eigenvalues, eigenvectors = eigsh(L_sparse, k=k, which=which)
            idx = np.argsort(eigenvalues) if which == "SM" else np.argsort(-eigenvalues)
            self.eigv = eigenvalues[idx]
            self.eigV = eigenvectors[:, idx]
        else:
            # Non-symmetric (rw): smallest/largest real part, real-cast & sorted
            from scipy.sparse.linalg import eigs
            which_ns = "SR" if which == "SM" else "LR"
            w, V = eigs(L_sparse, k=k, which=which_ns)
            self.eigv, self.eigV = _real_cast_sorted_gt(w, V)
            if which != "SM":
                self.eigv = self.eigv[::-1]
                self.eigV = self.eigV[:, ::-1]
        self._spectrum_laplacian_type = laplacian_type
        return self.eigv, self.eigV

    # ------------------------------------------------------------------
    # Adjacency spectrum methods
    # ------------------------------------------------------------------

    def compute_adjacency_spectrum_weigV(
        self,
        backend: Literal["numpy", "scipy", "cupy"] = "numpy",
        **kwargs: Any,
    ) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
        """Compute full eigendecomposition of the signed adjacency matrix.

        Parameters
        ----------
        backend : str
            ``'numpy'``, ``'scipy'``, or ``'cupy'``.
        **kwargs
            Ignored (NX compat: ``transpose``, ``flip_to_pos``, ``typf``).

        Returns
        -------
        tuple[ndarray, ndarray]
            ``(adj_eigv, adj_eigV)`` — eigenvalues sorted ascending,
            eigenvectors column-major.
        """
        A = self.get_signed_adjacency()

        if backend == "cupy":
            try:
                import cupy as cp
                eigv, eigV = cp.linalg.eigh(cp.asarray(A))
                eigv, eigV = cp.asnumpy(eigv), cp.asnumpy(eigV)
            except ImportError:
                eigv, eigV = np.linalg.eigh(A)
        elif backend == "scipy":
            from scipy.linalg import eigh
            eigv, eigV = eigh(A)
        else:
            eigv, eigV = np.linalg.eigh(A)

        idx = np.argsort(eigv)
        self.adj_eigv = eigv[idx]
        self.adj_eigV = eigV[:, idx]
        return self.adj_eigv, self.adj_eigV

    def compute_k_adj_eigvV(
        self,
        k: int = 1,
        which: str = "LM",
        **kwargs: Any,
    ) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
        """Compute k eigenvalues/eigenvectors of the adjacency matrix.

        Parameters
        ----------
        k : int
            Number of eigenpairs.
        which : str
            ``'LM'`` largest magnitude, ``'SM'`` smallest magnitude.
        **kwargs
            Ignored (NX compat).

        Returns
        -------
        tuple[ndarray, ndarray]
            ``(adj_eigv, adj_eigV)`` with k values.
        """
        A = self.get_signed_adjacency()
        A_sparse = csr_matrix(A)

        eigv, eigV = eigsh(A_sparse, k=k, which=which)
        idx = np.argsort(-eigv) if which == "LM" else np.argsort(eigv)
        self.adj_eigv = eigv[idx]
        self.adj_eigV = eigV[:, idx]
        return self.adj_eigv, self.adj_eigV

    # ------------------------------------------------------------------
    # Eigenvector format helpers
    # ------------------------------------------------------------------

    def make_eigV_transposed(self) -> None:
        """Transpose eigV to row-major (M x N) format."""
        if self.eigV is not None:
            self.eigV = self.eigV.T
            self._eigV_is_transposed = True

    def make_eigV_column_major(self) -> None:
        """Transpose eigV to column-major (N x M) format."""
        if self.eigV is not None:
            self.eigV = self.eigV.T
            self._eigV_is_transposed = False

    def make_adj_eigV_transposed(self) -> None:
        """Transpose adj_eigV to row-major."""
        if hasattr(self, "adj_eigV") and self.adj_eigV is not None:
            self.adj_eigV = self.adj_eigV.T
            self._adj_eigV_is_transposed = True

    def make_adj_eigV_column_major(self) -> None:
        """Transpose adj_eigV to column-major."""
        if hasattr(self, "adj_eigV") and self.adj_eigV is not None:
            self.adj_eigV = self.adj_eigV.T
            self._adj_eigV_is_transposed = False

    def make_rescaled_signed_laplacian(self, MODE: str = "field") -> None:
        """Compute rescaled signed Laplacian: L - lambda_0 * I.

        Parameters
        ----------
        MODE : str
            ``'field'`` subtracts smallest eigenvalue once.
            ``'double'`` subtracts twice (re-computes residual).
        """
        if self.eigv is None:
            self.compute_laplacian_spectrum()

        L = self.get_signed_laplacian()
        if MODE == "field":
            self.resLp = L - self.eigv[0] * np.eye(self.N)
        elif MODE == "double":
            from scipy.linalg import eigvalsh
            self.resLp = L - self.eigv[0] * np.eye(self.N)
            new_eigv0 = eigvalsh(self.resLp, subset_by_index=[0, 0])
            self.resLp = self.resLp - new_eigv0 * np.eye(self.N)

    def get_eigV(
        self,
        which: int = 0,
        binarize: bool = False,
    ) -> NDArray[np.floating]:
        """
        Get a specific eigenvector.

        Parameters
        ----------
        which : int, default 0
            Index of eigenvector (0 = smallest eigenvalue).
        binarize : bool, default False
            If True, return sign(eigenvector).

        Returns
        -------
        ndarray
            The eigenvector (or its binarization).

        Raises
        ------
        ValueError
            If eigenvectors haven't been computed.
        """
        if self.eigV is None:
            raise ValueError(
                "Eigenvectors not computed. "
                "Call compute_laplacian_spectrum_weigV() first."
            )

        eigvec = self.eigV[:, which]

        if binarize:
            return np.sign(eigvec).astype(np.int8)
        return eigvec

    def get_eigV_binarized(self, which: int = 0) -> NDArray[np.int8]:
        """
        Get binarized eigenvector (sign pattern).

        Parameters
        ----------
        which : int, default 0
            Index of eigenvector.

        Returns
        -------
        ndarray
            Array of +1/-1 values.
        """
        return self.get_eigV(which, binarize=True)

    def get_signed_laplacian_embedding(
        self,
        k: int = 2,
    ) -> NDArray[np.floating]:
        """
        Get k-dimensional spectral embedding from signed Laplacian.

        Parameters
        ----------
        k : int, default 2
            Number of dimensions (eigenvectors to use).

        Returns
        -------
        ndarray
            N x k embedding matrix.
        """
        if self.eigV is None or self.eigV.shape[1] < k:
            self.compute_laplacian_spectrum_weigV()

        # Use first k non-trivial eigenvectors (skip index 0 if it's trivial)
        return self.eigV[:, :k]

    # ------------------------------------------------------------------
    # Dynamics-critical methods (engine-agnostic interface)
    # ------------------------------------------------------------------

    def _build_neighbor_cache(self) -> None:
        """Pre-compute neighbor lists for all nodes (one-time GT traversal)."""
        sign_prop = self.G.edge_properties["sign"]
        cache: dict[int, list[tuple[int, float]]] = {}
        for node in range(self.N):
            v = self.G.vertex(node)
            neighbors: list[tuple[int, float]] = []
            for e in v.out_edges():
                nn = int(e.target()) if int(e.source()) == node else int(e.source())
                neighbors.append((nn, float(sign_prop[e])))
            cache[node] = neighbors
        self._neighbor_cache = cache

    def invalidate_neighbor_cache(self) -> None:
        """Invalidate cached neighbors (call after edge sign changes)."""
        self._neighbor_cache = None

    def get_neighbors_with_weights(self, node: int) -> list[tuple[int, float]]:
        """Return neighbors of *node* with signed edge weights.

        Uses a lazily-initialized cache to avoid repeated GT C++ traversals
        during dynamics simulations.

        Parameters
        ----------
        node : int
            Node index.

        Returns
        -------
        list[tuple[int, float]]
            List of ``(neighbor_index, sign_weight)`` pairs.
        """
        if self._neighbor_cache is None:
            self._build_neighbor_cache()
        return self._neighbor_cache[node]

    def get_edges_with_weights(self) -> list[tuple[int, int, float]]:
        """Return all edges with their signed weights.

        Returns
        -------
        list[tuple[int, int, float]]
            List of ``(u, v, weight)`` triples.
        """
        sign_prop = self.G.edge_properties["sign"]
        return [
            (int(e.source()), int(e.target()), float(sign_prop[e]))
            for e in self.G.edges()
        ]

    def get_edge_arrays(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (sources, targets, weights) as numpy arrays.

        Uses GT's bulk C++ accessors (``G.get_edges()`` and
        ``ep["sign"].a``) so that no per-edge Python loop is needed.

        Returns
        -------
        sources : ndarray[int64]
            Source vertex indices, shape ``(E,)``.
        targets : ndarray[int64]
            Target vertex indices, shape ``(E,)``.
        weights : ndarray[float64]
            Signed edge weights, shape ``(E,)``.
        """
        edges = self.G.get_edges()  # (E, 2) int64
        weights = self.G.ep["sign"].a.astype(np.float64)  # (E,)
        return edges[:, 0].copy(), edges[:, 1].copy(), weights

    def get_neighbor_indices(self, node: int) -> list[int]:
        """Return neighbor indices without weights.

        Parameters
        ----------
        node : int
            Node index.

        Returns
        -------
        list[int]
            Neighbor node indices.
        """
        v = self.G.vertex(node)
        return [int(n) for n in v.out_neighbors()]

    def invalidate_cache(self) -> None:
        """Clear cached matrices, spectral data, and neighbor cache."""
        self._slp = None
        self._adj = None
        self._eigv = None
        self._eigV = None
        self._spectrum_laplacian_type = None
        self._neighbor_cache = None

    def get_frustration_index(self) -> float:
        """
        Compute frustration index (fraction of frustrated triangles).

        Returns
        -------
        float
            Fraction of triangles that are frustrated (odd # negative edges).
        """
        # Count triangles and frustrated triangles
        # A triangle is frustrated if it has 1 or 3 negative edges
        sign_prop = self.G.edge_properties["sign"]
        n_triangles = 0
        n_frustrated = 0

        for v in self.G.vertices():
            neighbors = list(v.out_neighbors())
            for i, n1 in enumerate(neighbors):
                for n2 in neighbors[i + 1:]:
                    # Check if n1-n2 edge exists (completing triangle)
                    e12 = self.G.edge(n1, n2)
                    if e12 is not None:
                        n_triangles += 1
                        # Count negative edges in triangle
                        e_v_n1 = self.G.edge(v, n1)
                        e_v_n2 = self.G.edge(v, n2)
                        neg_count = sum([
                            sign_prop[e_v_n1] == -1,
                            sign_prop[e_v_n2] == -1,
                            sign_prop[e12] == -1,
                        ])
                        if neg_count % 2 == 1:  # 1 or 3 negative edges
                            n_frustrated += 1

        # Each triangle counted 3 times (once per vertex)
        n_triangles //= 3
        n_frustrated //= 3

        return n_frustrated / n_triangles if n_triangles > 0 else 0.0

    # ------------------------------------------------------------------
    # Eigenvector accessors (NX-compatible API)
    # ------------------------------------------------------------------

    def _get_eigvec(self, which: int = 0) -> NDArray[np.floating]:
        """Return eigenvector ``which`` from column-major eigV."""
        if self.eigV is None:
            raise ValueError(
                "Eigenvectors not computed. "
                "Call compute_k_eigvV() or compute_laplacian_spectrum_weigV() first."
            )
        if self.eigV.ndim == 1:
            if which != 0:
                raise IndexError(f"Only 1 eigenvector available, got which={which}")
            return self.eigV
        if self.eigV.shape[1] <= which:
            raise IndexError(
                f"Only {self.eigV.shape[1]} eigenvectors computed, "
                f"but index {which} requested."
            )
        return self.eigV[:, which]

    def get_eigV_check(
        self,
        which: int = 0,
        binarize: bool = False,
        **kwargs: Any,
    ) -> NDArray:
        """Get eigenvector ``which``, computing if not yet available.

        Parameters
        ----------
        which : int
            Eigenvector index (0 = Fiedler vector).
        binarize : bool
            If True, return sign(eigvec) in {-1, +1}.
        **kwargs
            Forwarded to ``compute_k_eigvV`` if recomputation is needed.
        """
        need_compute = (
            self.eigV is None
            or (self.eigV.ndim == 2 and self.eigV.shape[1] <= which)
        )
        if need_compute:
            self.compute_k_eigvV(k=which + 1)

        vec = self._get_eigvec(which)
        if binarize:
            return np.sign(vec).astype(np.int8)
        return vec

    def get_eigV_bin_check(
        self,
        which: int = 0,
        reshaped: bool = False,
        **kwargs: Any,
    ) -> NDArray[np.int8]:
        """Get binarized eigenvector with automatic computation.

        Convenience wrapper: ``get_eigV_check(which, binarize=True)``.
        NX-compatible alias used by ``BinDynSys.init_s()`` for
        ``init_cond='ground_state_*'``.

        Parameters
        ----------
        which : int
            Eigenvector index (0 = smallest eigenvalue).
        reshaped : bool
            Ignored (GT has no ``syshape``).  Accepted for NX API compat.
        **kwargs
            Forwarded to :meth:`get_eigV_check`.
        """
        return self.get_eigV_check(which, binarize=True, **kwargs)

    def get_eigV_bin_check_list(
        self,
        custom_slice: slice = slice(None),
        reshaped: bool = False,
        asarray: bool = False,
        **kwargs: Any,
    ) -> Union[List[NDArray], NDArray]:
        """Get multiple binarized eigenvectors.

        Parameters
        ----------
        custom_slice : slice
            ``slice.stop`` sets the number of eigenvectors.
            ``slice(None)`` returns all N.
        reshaped : bool
            Ignored (GT has no ``syshape``).
        asarray : bool
            If True, return stacked (k, N) array.
        **kwargs
            Forwarded to :meth:`get_eigV_check`.

        Returns
        -------
        list[ndarray] or ndarray
            Binarized eigenvectors in {-1, +1}.
        """
        stop = custom_slice.stop if custom_slice.stop is not None else self.N
        vecs = [self.get_eigV_check(i, binarize=True, **kwargs) for i in range(stop)]
        if asarray:
            return np.array(vecs)
        return vecs

    def load_eigV_on_graph(
        self,
        which: int = 0,
        on_g: Optional[str] = None,
        binarize: bool = False,
    ) -> None:
        """Ensure eigenvector is computed (GT does not store on node attrs).

        NX stores eigenvectors as node attributes on the graph object;
        GT keeps them in ``self.eigV`` and doesn't need that step.
        This method just ensures the eigenvector has been computed.

        Parameters
        ----------
        which : int
            Eigenvector index.
        on_g : str, optional
            Ignored (GT single representation).
        binarize : bool
            If True, also binarize (for validation; result is discarded).
        """
        self.get_eigV_check(which, binarize=binarize)

    def compute_laplacian_spectrum(
        self,
        typf: type = np.float64,
        backend: Optional[str] = None,
        keep_sparse: Optional[bool] = None,
        verbose: bool = False,
        laplacian_type: str = SG_LAPL_DEFAULT_TYPE,
    ) -> None:
        """Compute eigenvalues only (NX-compatible signature).

        Wraps :meth:`get_laplacian_spectrum` and caches the result
        in ``self.eigv``.  This method exists so that
        ``_ensure_spectrum()`` in the shared info-theory mixin
        can call ``self.compute_laplacian_spectrum(...)`` uniformly.

        Parameters
        ----------
        typf : type
            Ignored (GT always uses float64).
        backend : str, optional
            ``'numpy'``, ``'scipy'``, or ``'cupy'``.
        keep_sparse : bool, optional
            Ignored (GT always uses dense).
        verbose : bool
            Ignored.
        laplacian_type : str, default 'signed'
            Which Laplacian to diagonalize ('signed', 'sym', 'rw').
        """
        be = backend or "numpy"
        cached = getattr(self, "_eigv", None)
        if cached is not None and getattr(self, "_spectrum_laplacian_type", None) == laplacian_type:
            return
        self.eigv = self.get_laplacian_spectrum(backend=be, laplacian_type=laplacian_type)
        self._spectrum_laplacian_type = laplacian_type

    def get_sgspect_basis(
        self,
        max_factor: int = 2,
        start: int = 0,
        step: int = 1,
        normalized: bool = False,
    ) -> NDArray:
        """Return spectral basis matrix from eigenvectors.

        Parameters
        ----------
        max_factor : int
            Fraction of N to use as upper limit: ``N // max_factor``.
        start : int
            First eigenvector index.
        step : int
            Step between eigenvector indices.
        normalized : bool
            If True, L2-normalize each basis row.

        Returns
        -------
        ndarray, shape (k, N)
            Row-major basis matrix.
        """
        stop = int(self.N // max_factor)
        indices = range(start, stop, step)
        if not indices:
            return np.empty((0, self.N))
        # Ensure enough eigenvectors are computed
        max_idx = max(indices)
        if self.eigV is None or (self.eigV.ndim == 2 and self.eigV.shape[1] <= max_idx):
            self.compute_k_eigvV(k=max_idx + 1)

        basis = np.array([self._get_eigvec(i) for i in indices])
        if normalized:
            norms = np.linalg.norm(basis, axis=1, keepdims=True)
            norms = np.where(norms == 0, 1.0, norms)
            basis = basis / norms
        return basis

    # ------------------------------------------------------------------
    # Eigenvector clustering (scipy-based, no NX dependency)
    # ------------------------------------------------------------------

    def _find_components_in_subgraph(
        self, node_mask: NDArray[np.bool_]
    ) -> list[set[int]]:
        """Find connected components among nodes where *node_mask* is True.

        Uses the full signed adjacency matrix and scipy sparse connected
        components — no NX graph objects involved.

        Parameters
        ----------
        node_mask : ndarray of bool, shape (N,)
            True for nodes to include in the subgraph.

        Returns
        -------
        list[set[int]]
            Connected components, sorted by size descending.
        """
        indices = np.where(node_mask)[0]
        if len(indices) == 0:
            return []

        # Build adjacency sub-matrix for selected nodes
        adj_full = self.get_signed_adjacency()  # dense (N, N)
        # We only care about connectivity (ignore signs)
        sub_adj = np.abs(adj_full[np.ix_(indices, indices)])
        sub_sparse = csr_matrix(sub_adj)

        n_comp, labels = sp_connected_components(sub_sparse, directed=False)

        clusters: list[set[int]] = [set() for _ in range(n_comp)]
        for local_idx, label in enumerate(labels):
            clusters[label].add(int(indices[local_idx]))

        # Sort by size descending, drop empties
        clusters = sorted(
            (c for c in clusters if c), key=len, reverse=True
        )
        return clusters

    def make_clustersYN(
        self,
        which: int = 0,
        binarize: bool = True,
        val: int = 1,
    ) -> None:
        """Partition nodes by eigenvector sign and find connected components.

        Parameters
        ----------
        which : int
            Eigenvector index to partition on.
        binarize : bool
            If True, use sign of eigenvector.
        val : int
            Value considered "Y" partition (default +1).

        Sets
        ----
        self.clustersY : list[set[int]]
            Components where eigenvec matches *val*, sorted by size.
        self.clustersN : list[set[int]]
            Components where eigenvec does NOT match *val*.
        """
        eigvec = self.get_eigV_check(which, binarize=binarize)

        mask_y = eigvec == val
        mask_n = ~mask_y

        self.clustersY = self._find_components_in_subgraph(mask_y)
        self.clustersN = self._find_components_in_subgraph(mask_n)
        self.numClustersY = len(self.clustersY)
        self.numClustersN = len(self.clustersN)

        if self.clustersY:
            self.biggestClSet = self.clustersY
            self.numClustersBig = len(self.clustersY)
            self.gc = self.clustersY[0]
        else:
            self.biggestClSet = []
            self.numClustersBig = 0
            self.gc = None

    make_eigVclustersYN = make_clustersYN  # alias for NX API compat

    def get_eigV_cluster_sizes(
        self,
        which: int = 0,
        binarize: bool = True,
        val: int = 1,
    ) -> list[int]:
        """Return cluster sizes from eigenvector partitioning.

        Parameters
        ----------
        which : int
            Eigenvector index.
        binarize : bool
            Use sign of eigenvector.
        val : int
            Value for the "Y" partition.

        Returns
        -------
        list[int]
            Sorted cluster sizes (descending).
        """
        # Recompute clusters each time (stateless, avoids cache staleness)
        self.make_clustersYN(which=which, binarize=binarize, val=val)
        return sorted((len(c) for c in self.clustersY), reverse=True)

    def get_cluster_distribution(
        self,
        which: int = 0,
        binarize: bool = True,
        val: int = 1,
        **kwargs: Any,
    ) -> dict[int, int]:
        """Return cluster size distribution as ``{size: count}``.

        Parameters
        ----------
        which : int
            Eigenvector index.
        binarize : bool
            Use sign of eigenvector.
        val : int
            Value for "Y" partition.

        Returns
        -------
        dict[int, int]
            Mapping from cluster size to frequency.
        """
        cl_sizes = self.get_eigV_cluster_sizes(which, binarize, val)
        return {size: cl_sizes.count(size) for size in set(cl_sizes)}

    def compute_pinf(
        self,
        which: int = 0,
        val: int = 1,
        **kwargs: Any,
    ) -> None:
        """Compute infinite-cluster probability (Pinf).

        Parameters
        ----------
        which : int
            Eigenvector index.
        val : int
            Partition value.

        Sets
        ----
        self.Pinf : float
            Largest cluster size / N.
        self.Pinf_var : float
            Pinf variance estimate.
        """
        cl_sizes = np.array(self.get_eigV_cluster_sizes(which, True, val))
        if len(cl_sizes) == 0:
            self.Pinf = 0.0
            self.Pinf_var = 0.0
        else:
            self.Pinf = float(cl_sizes[0]) / self.N
            denominator = np.sum(cl_sizes) - cl_sizes[0]
            if denominator > 0:
                self.Pinf_var = float(
                    np.sum(cl_sizes @ cl_sizes - cl_sizes[0] ** 2)
                    / denominator
                )
            else:
                self.Pinf_var = 0.0

        if not hasattr(self, "Pinf_dict"):
            self.Pinf_dict: dict[int, tuple[float, float]] = {}
        self.Pinf_dict[which] = (self.Pinf, self.Pinf_var)

    def handle_no_clust(self, NoClust: int) -> int | None:
        """Check cluster availability (NX API compat).

        Parameters
        ----------
        NoClust : int
            Number of clusters requested.

        Returns
        -------
        int or None
            Validated cluster count, or None if no clusters exist.
        """
        num_big = getattr(self, "numClustersBig", 0)
        if num_big == 0:
            return None
        if NoClust > num_big:
            return num_big
        return NoClust

    # ------------------------------------------------------------------
    # Mixin imports — engine-agnostic methods from graphs._shared
    # ------------------------------------------------------------------

    # Info theory (Shannon & Renyi entropy, specific heat)
    from ._infotheory import compute_signed_laplacian_entropy
    from ._infotheory import compute_renyi_entropy_profile
    from ._infotheory import get_entropy, get_specific_heat
    from ._infotheory import get_entropy_derivative  # DEPRECATED
    from ._infotheory import get_renyi_results

    # Dynamics energy (RBIM, spherical SK from eigenvectors)
    from ._dynamics import compute_rbim_energy_eigV
    from ._dynamics import compute_rbim_energy_eigV_all
    from ._dynamics import get_rbim_energy_eigV
    from ._dynamics import get_all_rbim_energy_eigV
    from ._dynamics import compute_sksph_energy_eigV
    from ._dynamics import compute_sksph_energy_eigV_all
    from ._dynamics import get_sksph_energy_eigV
    from ._dynamics import get_all_sksph_energy_eigV

    # Order parameters (spectral gap)
    from ._ordparams import compute_gap
    from ._ordparams import compute_gap_between
    from ._ordparams import get_gap

    # Quantum propagator
    from ._quantum import compute_quantum_propagator
    from ._quantum import quantum_walk_probabilities
    from ._quantum import quantum_observables_time_series

    # Graph operations (NX _ongraph.py parity)
    from ._ongraph import check_Ne_flips
    from ._ongraph import set_edges_random_normal
    from ._ongraph import load_vec_on_nodes
    from ._ongraph import set_node_attributes
    from ._ongraph import get_random_edges_from_set

    # Partitioning (NX _partitioning.py parity)
    from ._partitioning import get_subgraph_from_nodes
    from ._partitioning import get_nodes_subgraph_by_kv
    from ._partitioning import make_graphYN
    from ._partitioning import make_connected_component_by_edge
    from ._partitioning import get_ferroAntiferro_regions

    # File I/O — exports
    from ._exports import export_eigV_all
    from ._exports import export_adj_bin
    from ._exports import export_ising_clust

    # File I/O — loaders
    from ._loaders import load_eigV_all
    from ._loaders import set_edgel_from_bin

    # File I/O — cleaners
    from ._cleaners import remove_ising_clust_files
    from ._cleaners import remove_edgl_file
    from ._cleaners import remove_eigV_file
    from ._cleaners import remove_adj_file
    from ._cleaners import remove_exported_files
    from ._cleaners import clean_gclutil

    # Topology helpers (NX _topology.py parity)
    from ._topology import nodes_in
    from ._topology import get_node_attributes
    from ._topology import get_edge_data
    from ._topology import get_edge_color
    from ._topology import get_graph_neighbors
    from ._topology import get_adjacency_matrix_for
    from ._topology import get_degree_matrix_for
    from ._topology import get_laplacian_matrix_for
    from ._topology import get_signed_degree_matrix_for
    from ._topology import get_signed_laplacian_matrix_for

    # Representation stubs (NX multi-repr compat)
    from ._representations import upd_graph_matrices
    from ._representations import upd_edge_sets
    from ._representations import upd_Degree
    from ._representations import upd_GraphRepr_All
    from ._representations import upd_NodeMap
    from ._representations import upd_EdgeMap
    from ._representations import upd_ReprMaps
    from ._representations import zip_reprNodes
    from ._representations import zip_reprEdges

    @property
    def gcl(self):
        """Graph clustering utility (lazy-initialized NestedDict)."""
        if not hasattr(self, "_gcl"):
            from ....utils.tools import NestedDict
            self._gcl = NestedDict()
        return self._gcl

    @property
    def graph_clustering_utility(self):
        """Alias for gcl (NX compatibility)."""
        return self.gcl

    # ==================================================================
    # DynamicsGraphProtocol — path infrastructure & C backend I/O
    # ==================================================================

    def _ensure_paths(self) -> None:
        """Lazily initialise directory paths for C backend file I/O.

        Mirrors the NX path tree:
        ``path_sgdata / {graph,lrgsg,phtra,spect,ising,...} / syshapePth``

        Called on first access to any path property.  Subclasses that
        set ``self.syshapePth`` in their own ``__init__`` (e.g.
        ``Lattice3DGT``) will see that value used here.
        """
        if self._paths_initialised:
            return
        self._paths_initialised = True

        self.path_data = self._path_data or PATHDATA
        self._path_sgdata = self.path_data / Path(self._sgpathn)

        # syshapePth may have been set by a subclass already
        if not hasattr(self, "_syshapePth") or self._syshapePth is None:
            self._syshapePth = f"N={self.N}"

        for p in PATHN_GRAPH_LIST + PATHN_DYNAMICS_LIST:
            pfname = Path(p, self._syshapePth)
            setattr(self, f"path_{p}", self._path_sgdata / pfname)

    @property
    def syshapePth(self) -> str:
        """Parameter string used in directory/file naming."""
        if hasattr(self, "_syshapePth") and self._syshapePth is not None:
            return self._syshapePth
        return f"N={self.N}"

    @syshapePth.setter
    def syshapePth(self, value: str) -> None:
        self._syshapePth = value

    @property
    def path_sgdata(self) -> Path:
        """Root data directory for this graph instance."""
        self._ensure_paths()
        return self._path_sgdata

    def get_p_fname(
        self,
        who: str,
        out_suffix: str = "",
        ext: str = BIN,
    ) -> str | Path:
        """Build a parametric filename.

        Parameters
        ----------
        who : str
            Filename prefix (e.g. ``'s'``, ``'h'``, ``'edgelist'``).
        out_suffix : str
            Additional suffix appended before the extension.
        ext : str
            File extension (default ``'.bin'``).

        Returns
        -------
        Path
            Constructed filename (relative, no directory).
        """
        return build_p_fname(who, self.pflip, out_suffix=out_suffix, ext=ext)

    def _export_edgel_bin(self, exName: str = "") -> None:
        """Write edge list to binary file for C backend.

        Binary format per edge: ``(uint64 i, uint64 j, float64 w_ij)``
        — 24 bytes, matching the NX implementation exactly.

        Parameters
        ----------
        exName : str
            Suffix appended to the filename (typically the run id).
        """
        self._ensure_paths()
        fname = build_p_fname("edgelist", self.pflip, out_suffix=exName, ext=BIN)
        self.path_exp_edgl = self.path_graph / fname
        self.path_exp_edgl.parent.mkdir(parents=True, exist_ok=True)

        sign_prop = self.G.edge_properties["sign"]
        dtype = [("i", np.uint64), ("j", np.uint64), ("w_ij", np.float64)]
        edge_array = np.array(
            [(int(e.source()), int(e.target()), float(sign_prop[e]))
             for e in self.G.edges()],
            dtype=dtype,
        )
        with open(self.path_exp_edgl, "wb") as f:
            edge_array.tofile(f)

    def remove_exported_files(self) -> None:
        """Remove all exported binary files created by C backend I/O."""
        for attr in ("path_exp_edgl",):
            path = getattr(self, attr, None)
            if path is not None and os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"SignedGraphGT(N={self.N}, edges={self.num_edges}, "
            f"negative={neg}, pflip={self.pflip})"
        )


__all__ = ["SignedGraphGT"]
