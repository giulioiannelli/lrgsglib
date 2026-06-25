"""Engine-free graph specification — the single return contract for generators.

A :class:`GraphSpec` is a pure topology description (node count, integer edge
list, optional positions) with **no backend object**. It is what every graph
generator should return, regardless of how it is implemented (pure Python,
vectorized NumPy, or a compiled pybind/C++ extension).

The point of the spec is decoupling: a generator is written *once* and emits a
spec; the spec is then *materialized* into whichever engine the caller wants
(NetworkX or graph-tool), and converting a built graph from one engine to the
other is just "extract the spec, materialize into the other engine". One
generator therefore serves both engines and both conversion directions for free.

Layering
--------
This module operates on the **raw backend graph objects** (``networkx.Graph`` /
``graph_tool.Graph``), not on the ``SignedGraph`` wrappers, so it stays the
lowest reusable layer. Signs are intentionally *not* part of the spec: the
``SignedGraph`` subclasses assign the sign edge-property after materialization
(matching the current GT build, which is all-positive at construction and
flipped afterwards).

Optional engine-native fast path
--------------------------------
The default contract is "generator returns a :class:`GraphSpec`". A generator
*may additionally* expose a ``gt_native`` builder that constructs a
``graph_tool.Graph`` directly, used only when ``engine == 'gt'`` and the
NumPy-array round-trip has been profiled as a bottleneck for that specific
geometry. :func:`build` performs that dispatch. Most generators never need it —
``graph-tool``'s ``add_edge_list`` is already a vectorized C++ ingest.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Protocol, Tuple, Union, runtime_checkable

import numpy as np

from .._engine import GraphEngine

__all__ = [
    "GraphSpec",
    "materialize",
    "materialize_nx",
    "materialize_gt",
    "spec_from_nx",
    "spec_from_gt",
    "spec_from_graph",
    "detect_engine",
    "convert",
    "build",
    "SpecGenerator",
]

# Vertex-property / node-attribute key under which positions are stored, matching
# the existing engines (GT ``vertex_properties["pos"]`` and NX node attr "pos").
_POS_KEY = "pos"


@dataclass
class GraphSpec:
    """Engine-free description of a graph's topology.

    Parameters
    ----------
    n : int
        Number of vertices, labelled ``0 .. n-1``.
    edges : np.ndarray
        Integer array of shape ``(M, 2)``; each row is an undirected edge
        ``(u, v)`` with ``u, v in [0, n)``. No self-loops or duplicates are
        assumed to have been removed by the generator.
    pos : np.ndarray, optional
        Float array of shape ``(n, d)`` (``d`` = 2 or 3) with one position per
        vertex, or ``None`` if the geometry carries no coordinates.
    directed : bool, default False
        Whether the edge list should be materialized as a directed graph.

    Notes
    -----
    This is a topology-only primitive — it deliberately carries no edge signs or
    weights; those are applied by the ``SignedGraph`` wrapper after
    materialization.
    """

    n: int
    edges: np.ndarray
    pos: Optional[np.ndarray] = None
    directed: bool = False

    def __post_init__(self) -> None:
        self.edges = np.asarray(self.edges, dtype=np.int64).reshape(-1, 2)
        if self.pos is not None:
            self.pos = np.asarray(self.pos, dtype=np.float64)
            if self.pos.shape[0] != self.n:
                raise ValueError(
                    f"pos has {self.pos.shape[0]} rows but n={self.n}"
                )
        if self.edges.size and int(self.edges.max()) >= self.n:
            raise ValueError(
                f"edge endpoint {int(self.edges.max())} >= n={self.n}"
            )

    @property
    def num_edges(self) -> int:
        """Number of edges (rows of :attr:`edges`)."""
        return int(self.edges.shape[0])

    @classmethod
    def from_genresult(
        cls, result: Tuple[int, Any, Optional[np.ndarray]], directed: bool = False
    ) -> "GraphSpec":
        """Build a spec from the legacy ``(n_nodes, edges, pos)`` generator tuple.

        Bridges the existing ``graphs/gt/*/_generators.py`` ``GenResult`` tuples
        onto the shared contract without rewriting them.
        """
        n, edges, pos = result
        return cls(n=int(n), edges=np.asarray(edges, dtype=np.int64),
                   pos=pos, directed=directed)


@runtime_checkable
class SpecGenerator(Protocol):
    """A generator is any callable returning a :class:`GraphSpec`.

    It may optionally carry a ``gt_native`` attribute (a callable with the same
    signature returning a ``graph_tool.Graph`` directly) for the engine-native
    fast path; :func:`build` uses it when present and ``engine == 'gt'``.
    """

    def __call__(self, *args: Any, **kwargs: Any) -> GraphSpec: ...


# --------------------------------------------------------------------------- #
# Materialization: spec -> raw backend graph
# --------------------------------------------------------------------------- #
def materialize_gt(spec: GraphSpec, with_positions: bool = True):
    """Materialize a spec into a raw ``graph_tool.Graph``.

    Uses ``add_edge_list`` (vectorized C++ ingest) and ``set_2d_array`` for
    positions, so the cost is dominated by graph-tool's native code rather than
    Python. The returned graph carries no sign property — the ``SignedGraphGT``
    wrapper assigns that.
    """
    import graph_tool.all as gt

    G = gt.Graph(directed=spec.directed)
    G.add_vertex(spec.n)
    if spec.num_edges:
        G.add_edge_list(spec.edges)
    if with_positions and spec.pos is not None:
        vp = G.new_vertex_property("vector<double>")
        vp.set_2d_array(spec.pos.T)
        G.vertex_properties[_POS_KEY] = vp
    return G


def materialize_nx(spec: GraphSpec, with_positions: bool = True):
    """Materialize a spec into a raw ``networkx.Graph`` with integer nodes.

    Builds the adjacency as a SciPy CSR matrix and uses
    ``from_scipy_sparse_array`` — this avoids both the Python-level
    ``add_edges_from`` loop and the ``convert_node_labels_to_integers`` relabel
    that the legacy NX lattice path pays. All ``n`` vertices are present even if
    some are isolated.
    """
    import networkx as nx
    from scipy.sparse import coo_matrix

    cls = nx.DiGraph if spec.directed else nx.Graph
    if spec.num_edges:
        u, v = spec.edges[:, 0], spec.edges[:, 1]
        if spec.directed:
            rows, cols = u, v
        else:
            rows = np.concatenate([u, v])
            cols = np.concatenate([v, u])
        data = np.ones(rows.shape[0], dtype=np.int8)
        adj = coo_matrix((data, (rows, cols)), shape=(spec.n, spec.n)).tocsr()
        G = nx.from_scipy_sparse_array(
            adj, create_using=cls, edge_attribute=None)
    else:
        G = cls()
        G.add_nodes_from(range(spec.n))
    if with_positions and spec.pos is not None:
        nx.set_node_attributes(
            G, {i: tuple(spec.pos[i]) for i in range(spec.n)}, _POS_KEY)
    return G


def materialize(
    spec: GraphSpec,
    engine: Union[str, GraphEngine],
    with_positions: bool = True,
):
    """Materialize a spec into the backend graph for ``engine`` (``'nx'``/``'gt'``)."""
    engine = GraphEngine(engine) if isinstance(engine, str) else engine
    if engine is GraphEngine.NETWORKX:
        return materialize_nx(spec, with_positions=with_positions)
    if engine is GraphEngine.GRAPHTOOL:
        return materialize_gt(spec, with_positions=with_positions)
    raise ValueError(f"unsupported engine for materialization: {engine}")


# --------------------------------------------------------------------------- #
# Extraction: raw backend graph -> spec
# --------------------------------------------------------------------------- #
def spec_from_gt(G) -> GraphSpec:
    """Extract a :class:`GraphSpec` from a raw ``graph_tool.Graph``.

    Uses ``get_edges`` (returns an ``ndarray``) for a vectorized read; positions
    are read from ``vertex_properties["pos"]`` if present.
    """
    edges = np.asarray(G.get_edges(), dtype=np.int64)[:, :2]
    pos = None
    if _POS_KEY in G.vertex_properties:
        vp = G.vertex_properties[_POS_KEY]
        d = len(vp[G.vertex(0)]) if G.num_vertices() else 2
        pos = vp.get_2d_array(range(d)).T
    return GraphSpec(n=int(G.num_vertices()), edges=edges, pos=pos,
                     directed=bool(G.is_directed()))


def spec_from_nx(G) -> GraphSpec:
    """Extract a :class:`GraphSpec` from a raw ``networkx.Graph``.

    Requires integer node labels ``0 .. n-1`` (the materialized convention). Uses
    ``to_scipy_sparse_array`` for a vectorized edge read; positions come from the
    "pos" node attribute if present.
    """
    import networkx as nx
    from scipy.sparse import triu

    n = G.number_of_nodes()
    nodes = list(G.nodes())
    if nodes and (min(nodes) != 0 or max(nodes) != n - 1):
        raise ValueError(
            "spec_from_nx requires integer node labels 0..n-1; "
            "relabel before extracting a spec")
    adj = nx.to_scipy_sparse_array(G, nodelist=range(n), format="coo")
    if not G.is_directed():
        adj = triu(adj, k=1).tocoo()
    edges = np.column_stack([adj.row, adj.col]).astype(np.int64)
    pos = None
    pos_attr = nx.get_node_attributes(G, _POS_KEY)
    if pos_attr:
        pos = np.array([pos_attr[i] for i in range(n)], dtype=np.float64)
    return GraphSpec(n=int(n), edges=edges, pos=pos,
                     directed=bool(G.is_directed()))


def detect_engine(G) -> GraphEngine:
    """Return the :class:`GraphEngine` of a raw backend graph object."""
    import networkx as nx
    if isinstance(G, nx.Graph):
        return GraphEngine.NETWORKX
    try:
        import graph_tool.all as gt
        if isinstance(G, gt.Graph):
            return GraphEngine.GRAPHTOOL
    except ImportError:
        pass
    raise TypeError(f"unrecognized graph object of type {type(G)!r}")


def spec_from_graph(G) -> GraphSpec:
    """Extract a spec from a raw backend graph, auto-detecting the engine."""
    engine = detect_engine(G)
    return (spec_from_nx(G) if engine is GraphEngine.NETWORKX
            else spec_from_gt(G))


# --------------------------------------------------------------------------- #
# Conversion and generator-driven build
# --------------------------------------------------------------------------- #
def convert(G, to: Union[str, GraphEngine], with_positions: bool = True):
    """Convert a raw backend graph to another engine via the spec.

    Equivalent to ``materialize(spec_from_graph(G), to)``: extract the topology
    once, rebuild it in the target engine. Replaces per-engine object-walking
    converters with a single vectorized path. Edge signs/weights are *not*
    carried (this is a topology conversion); re-apply them on the wrapper.
    """
    return materialize(spec_from_graph(G), to, with_positions=with_positions)


def build(
    generator: SpecGenerator,
    engine: Union[str, GraphEngine],
    *args: Any,
    with_positions: bool = True,
    **kwargs: Any,
):
    """Run a generator and materialize its spec into ``engine``.

    Implements the optional engine-native fast path: if ``engine == 'gt'`` and
    ``generator`` exposes a ``gt_native`` callable, that builder is used directly
    (bypassing the array round-trip); otherwise the generator's
    :class:`GraphSpec` is materialized normally.
    """
    engine = GraphEngine(engine) if isinstance(engine, str) else engine
    native = getattr(generator, "gt_native", None)
    if engine is GraphEngine.GRAPHTOOL and callable(native):
        return native(*args, **kwargs)
    spec = generator(*args, **kwargs)
    if not isinstance(spec, GraphSpec):
        spec = GraphSpec.from_genresult(spec)
    return materialize(spec, engine, with_positions=with_positions)
