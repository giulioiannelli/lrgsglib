"""
TemporalGraphNX: Time-varying signed network with snapshot-based representation.

Inherits from SignedGraphNX, providing full signed graph functionality
(spectral analysis, entropy, clustering) on any materialized snapshot.

Example
-------
>>> from lrgsglib.graphs.nx import TemporalGraphNX
>>> tg = TemporalGraphNX(n_nodes=100, pflip=0.1, seed=42)
>>> tg.add_edge(0, 1, t_start=0.5, sign=1)
>>> tg.add_edge(1, 2, t_start=1.0, sign=-1)
>>> snapshot = tg.get_snapshot(time=0.75)
>>> sg = tg.get_signed_snapshot(time=1.5)
"""

from typing import Optional, Union, Sequence, Iterator
from collections import defaultdict
import networkx as nx
import numpy as np

from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


class TemporalGraphNX(SignedGraphNX):
    """
    Time-varying signed network with support for snapshots and edge streams.

    Inherits from SignedGraphNX so that all signed graph functionality
    (spectral analysis, entropy, clustering, dynamics) is available.
    The base class operates on an initial empty graph with n_nodes;
    use ``get_signed_snapshot(time)`` to obtain a fully materialized
    SignedGraphNX at a specific time point.

    Parameters
    ----------
    n_nodes : int, optional
        Number of nodes (fixed throughout time). If None, nodes are
        added dynamically.
    n_snapshots : int, optional
        If provided, pre-allocate storage for discrete snapshots.
    time_unit : str, default "continuous"
        Time representation: "continuous" for real-valued times,
        "discrete" for integer time steps.
    directed : bool, default False
        Whether edges are directed.
    pflip : float, default 0.0
        Probability of sign flip (fraction of negative edges).
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    n_nodes : int or None
        Number of nodes.
    temporal_edges : list
        List of (u, v, t_start, t_end, attrs) tuples.
    times : list
        Sorted list of all time points where changes occur.
    time_range : tuple
        (min_time, max_time) tuple.

    Notes
    -----
    Two temporal representations are supported:

    1. **Edge stream**: Each edge has a time interval [t_start, t_end)
    2. **Snapshots**: Discrete graph states at specific times

    Edge signs can vary over time via ``flip_random_edges()``.

    The inherited SignedGraphNX operates on the initial (empty) graph.
    For spectral analysis at a specific time, use::

        sg = tg.get_signed_snapshot(time=5.0)
        sg.compute_laplacian_spectrum_weigV()

    Example
    -------
    >>> tg = TemporalGraphNX(n_nodes=50, pflip=0.1, seed=42)
    >>> tg.add_temporal_edge(0, 1, t_start=0, t_end=10, sign=1)
    >>> tg.add_temporal_edge(1, 2, t_start=5, t_end=15, sign=-1)
    >>> G = tg.get_snapshot(7)
    >>> G.number_of_edges()
    2
    """

    sgpathn = "temporal_graph"

    def __init__(
        self,
        n_nodes: Optional[int] = None,
        n_snapshots: Optional[int] = None,
        time_unit: str = "continuous",
        directed: bool = False,
        pflip: float = 0.0,
        seed: Optional[int] = None,
    ) -> None:
        # Create initial graph with n_nodes (no edges)
        G = nx.DiGraph() if directed else nx.Graph()
        if n_nodes is not None:
            G.add_nodes_from(range(n_nodes))

        # Initialize SignedGraphNX base class
        super().__init__(
            G,
            pflip=pflip,
            seed=seed,
            make_dir_tree=False,
        )

        # Temporal-specific attributes
        self.n_nodes = n_nodes
        self.n_snapshots = n_snapshots
        self.time_unit = time_unit
        self.directed = directed

        # Edge stream: list of (u, v, t_start, t_end, attrs)
        self._temporal_edges: list[tuple] = []

        # Index edges by time for efficient queries
        self._edge_starts: dict[float, list[int]] = defaultdict(list)
        self._edge_ends: dict[float, list[int]] = defaultdict(list)

        # Cached snapshots
        self._temporal_snapshots: dict[float, nx.Graph] = {}

        # Track all time points
        self._times: set[float] = set()

        # Node set for temporal tracking (can grow if n_nodes is None)
        self._temporal_nodes: set[int] = set(range(n_nodes)) if n_nodes else set()

        # RNG for sign flipping
        self._temporal_rng = np.random.default_rng(seed)

        # Track sign flip events: list of (time, u, v, new_sign)
        self._sign_flips: list[tuple[float, int, int, int]] = []

    # ------------------------------------------------------------------
    # Temporal properties
    # ------------------------------------------------------------------

    @property
    def times(self) -> list[float]:
        """Return sorted list of all time points."""
        return sorted(self._times)

    @property
    def time_range(self) -> tuple[float, float]:
        """Return (min_time, max_time) tuple."""
        if not self._times:
            return (0.0, 0.0)
        return (min(self._times), max(self._times))

    @property
    def num_edges_total(self) -> int:
        """Return total number of temporal edges."""
        return len(self._temporal_edges)

    @property
    def temporal_edges(self) -> list[tuple]:
        """Return the temporal edge stream."""
        return self._temporal_edges

    # ------------------------------------------------------------------
    # Edge / node management
    # ------------------------------------------------------------------

    def add_node(self, node: int) -> None:
        """Add a node to the temporal graph."""
        self._temporal_nodes.add(node)

    def add_temporal_edge(
        self,
        u: int,
        v: int,
        t_start: float = 0.0,
        t_end: Optional[float] = None,
        sign: int = 1,
        weight: float = 1.0,
        **attrs,
    ) -> int:
        """
        Add a temporal edge.

        Parameters
        ----------
        u, v : int
            Edge endpoints.
        t_start : float, default 0.0
            Time when edge appears.
        t_end : float, optional
            Time when edge disappears. If None, edge persists forever.
        sign : int, default 1
            Edge sign (+1 or -1).
        weight : float, default 1.0
            Edge weight.
        **attrs
            Additional edge attributes.

        Returns
        -------
        int
            Index of the added edge.
        """
        self._temporal_nodes.add(u)
        self._temporal_nodes.add(v)

        if t_end is None:
            t_end = float("inf")

        edge_idx = len(self._temporal_edges)
        edge_attrs = {"sign": sign, "weight": weight, **attrs}
        self._temporal_edges.append((u, v, t_start, t_end, edge_attrs))

        self._edge_starts[t_start].append(edge_idx)
        self._edge_ends[t_end].append(edge_idx)

        self._times.add(t_start)
        if t_end != float("inf"):
            self._times.add(t_end)

        self._invalidate_temporal_cache(t_start, t_end)
        return edge_idx

    # Keep add_edge as alias for backward compatibility
    def add_edge(
        self,
        u: int,
        v: int,
        t_start: float = 0.0,
        t_end: Optional[float] = None,
        sign: int = 1,
        weight: float = 1.0,
        **attrs,
    ) -> int:
        """Alias for ``add_temporal_edge()`` (backward compatibility)."""
        return self.add_temporal_edge(
            u, v, t_start=t_start, t_end=t_end,
            sign=sign, weight=weight, **attrs,
        )

    def remove_edge(self, u: int, v: int, time: float) -> bool:
        """
        Remove an edge at a specific time (sets t_end to given time).

        Parameters
        ----------
        u, v : int
            Edge endpoints.
        time : float
            Time at which to remove the edge.

        Returns
        -------
        bool
            True if edge was found and removed.
        """
        for idx, (eu, ev, t_start, t_end, attrs) in enumerate(self._temporal_edges):
            if (eu == u and ev == v) or (not self.directed and eu == v and ev == u):
                if t_start <= time < t_end:
                    self._temporal_edges[idx] = (eu, ev, t_start, time, attrs)
                    self._times.add(time)
                    self._invalidate_temporal_cache(time, t_end)
                    return True
        return False

    # ------------------------------------------------------------------
    # Snapshots
    # ------------------------------------------------------------------

    def add_snapshot(
        self,
        time: float,
        edges: Optional[Sequence[tuple]] = None,
        graph: Optional[nx.Graph] = None,
    ) -> None:
        """
        Add a complete snapshot at a specific time.

        Parameters
        ----------
        time : float
            Snapshot time.
        edges : sequence of tuples, optional
            List of (u, v) or (u, v, attrs) tuples.
        graph : nx.Graph, optional
            Complete graph snapshot.
        """
        if graph is not None:
            self._temporal_snapshots[time] = graph.copy()
            self._temporal_nodes.update(graph.nodes())
        elif edges is not None:
            G = nx.DiGraph() if self.directed else nx.Graph()
            for edge in edges:
                if len(edge) == 2:
                    G.add_edge(edge[0], edge[1], sign=1, weight=1.0)
                elif len(edge) == 3:
                    G.add_edge(edge[0], edge[1], **edge[2])
                else:
                    G.add_edge(edge[0], edge[1])
            self._temporal_snapshots[time] = G
            self._temporal_nodes.update(G.nodes())

        self._times.add(time)

    def get_snapshot(self, time: float) -> nx.Graph:
        """
        Get a static graph snapshot at the given time.

        Parameters
        ----------
        time : float
            Time at which to extract snapshot.

        Returns
        -------
        nx.Graph or nx.DiGraph
            Static graph with edges active at the given time.
        """
        if time in self._temporal_snapshots:
            return self._temporal_snapshots[time].copy()

        G = nx.DiGraph() if self.directed else nx.Graph()
        G.add_nodes_from(self._temporal_nodes)

        for u, v, t_start, t_end, attrs in self._temporal_edges:
            if t_start <= time < t_end:
                G.add_edge(u, v, **attrs)

        return G

    def get_signed_snapshot(self, time: float) -> SignedGraphNX:
        """
        Get a full SignedGraphNX snapshot at the given time.

        The returned object has all SignedGraphNX methods available
        (spectral analysis, entropy, clustering, etc.).

        Parameters
        ----------
        time : float
            Time at which to extract snapshot.

        Returns
        -------
        SignedGraphNX
            Signed graph with edges and signs at the given time.
        """
        G = self.get_snapshot(time)

        # Apply sign flips up to this time
        for flip_time, u, v, new_sign in self._sign_flips:
            if flip_time <= time and G.has_edge(u, v):
                G[u][v]["sign"] = new_sign
                G[u][v]["weight"] = new_sign

        return SignedGraphNX(G, pflip=self.pflip, seed=self.seed,
                             make_dir_tree=False)

    def get_aggregated_graph(
        self,
        t_start: Optional[float] = None,
        t_end: Optional[float] = None,
        aggregation: str = "union",
    ) -> nx.Graph:
        """
        Get aggregated graph over a time window.

        Parameters
        ----------
        t_start : float, optional
            Start of time window (default: min time).
        t_end : float, optional
            End of time window (default: max time).
        aggregation : str, default "union"
            Aggregation method:
            - "union": Include edge if active at any time
            - "intersection": Include edge only if always active
            - "weighted": Weight by duration active

        Returns
        -------
        nx.Graph
            Aggregated graph.
        """
        if t_start is None:
            t_start = self.time_range[0]
        if t_end is None:
            t_end = self.time_range[1]

        G = nx.DiGraph() if self.directed else nx.Graph()
        G.add_nodes_from(self._temporal_nodes)

        if aggregation == "union":
            for u, v, ts, te, attrs in self._temporal_edges:
                if ts < t_end and te > t_start:
                    if G.has_edge(u, v):
                        G[u][v]["weight"] += attrs.get("weight", 1.0)
                    else:
                        G.add_edge(u, v, **attrs)

        elif aggregation == "intersection":
            edge_counts = defaultdict(int)
            total_times = len(self._times)
            for t in self._times:
                if t_start <= t < t_end:
                    for u, v, ts, te, attrs in self._temporal_edges:
                        if ts <= t < te:
                            edge_counts[(u, v)] += 1
            for (u, v), count in edge_counts.items():
                if count == total_times:
                    for eu, ev, ts, te, attrs in self._temporal_edges:
                        if (eu, ev) == (u, v):
                            G.add_edge(u, v, **attrs)
                            break

        elif aggregation == "weighted":
            duration = t_end - t_start
            for u, v, ts, te, attrs in self._temporal_edges:
                overlap_start = max(ts, t_start)
                overlap_end = min(te, t_end)
                if overlap_start < overlap_end:
                    weight = (overlap_end - overlap_start) / duration
                    if G.has_edge(u, v):
                        G[u][v]["weight"] += weight
                    else:
                        G.add_edge(u, v, weight=weight,
                                   sign=attrs.get("sign", 1))

        return G

    # ------------------------------------------------------------------
    # Iteration
    # ------------------------------------------------------------------

    def iter_snapshots(
        self,
        step: Optional[float] = None,
    ) -> Iterator[tuple[float, nx.Graph]]:
        """
        Iterate over snapshots.

        Parameters
        ----------
        step : float, optional
            Time step between snapshots. If None, yields at each
            time point where a change occurs.

        Yields
        ------
        tuple[float, nx.Graph]
            (time, snapshot) pairs.
        """
        if step is None:
            for t in self.times:
                yield t, self.get_snapshot(t)
        else:
            t_min, t_max = self.time_range
            t = t_min
            while t <= t_max:
                yield t, self.get_snapshot(t)
                t += step

    # ------------------------------------------------------------------
    # Activity queries
    # ------------------------------------------------------------------

    def get_edge_activity(self, u: int, v: int) -> list[tuple[float, float]]:
        """
        Get time intervals when an edge is active.

        Parameters
        ----------
        u, v : int
            Edge endpoints.

        Returns
        -------
        list of (t_start, t_end) tuples
            Time intervals when edge is active.
        """
        intervals = []
        for eu, ev, ts, te, _ in self._temporal_edges:
            if (eu == u and ev == v) or (not self.directed and eu == v and ev == u):
                intervals.append((ts, te))
        return sorted(intervals)

    def get_node_activity(self, node: int) -> list[tuple[float, float]]:
        """
        Get time intervals when a node has any edges.

        Parameters
        ----------
        node : int
            Node ID.

        Returns
        -------
        list of (t_start, t_end) tuples
            Time intervals when node has at least one edge.
        """
        intervals = []
        for u, v, ts, te, _ in self._temporal_edges:
            if u == node or v == node:
                intervals.append((ts, te))
        return self._merge_intervals(intervals)

    # ------------------------------------------------------------------
    # Sign flipping (merged from TemporalSignedGraphNX)
    # ------------------------------------------------------------------

    def flip_random_edges(
        self,
        time: float,
        fraction: Optional[float] = None,
        edges: Optional[Sequence[tuple]] = None,
    ) -> list[tuple]:
        """
        Flip signs of random edges at a specific time.

        Parameters
        ----------
        time : float
            Time at which to flip edges.
        fraction : float, optional
            Fraction of edges to flip (default: self.pflip).
        edges : sequence of (u, v) tuples, optional
            Specific edges to flip. If None, select randomly.

        Returns
        -------
        list of (u, v) tuples
            Edges that were flipped.
        """
        if fraction is None:
            fraction = self.pflip

        G = self.get_snapshot(time)

        if edges is None:
            all_edges = list(G.edges())
            if not all_edges:
                return []
            n_flip = max(1, int(len(all_edges) * fraction))
            flip_indices = self._temporal_rng.choice(
                len(all_edges), size=min(n_flip, len(all_edges)), replace=False
            )
            edges = [all_edges[i] for i in flip_indices]

        flipped = []
        for u, v in edges:
            if G.has_edge(u, v):
                old_sign = G[u][v].get("sign", 1)
                new_sign = -old_sign
                self._sign_flips.append((time, u, v, new_sign))
                flipped.append((u, v))

        self._times.add(time)
        self._invalidate_temporal_cache(time, float("inf"))
        return flipped

    def get_frustration_evolution(
        self,
        times: Optional[Sequence[float]] = None,
    ) -> list[tuple[float, float]]:
        """
        Compute frustration index evolution over time.

        Parameters
        ----------
        times : sequence of float, optional
            Times at which to compute frustration. Default: all event times.

        Returns
        -------
        list of (time, frustration) tuples
            Frustration index at each time point.
        """
        if times is None:
            times = self.times

        frustrations = []
        for t in times:
            sg = self.get_signed_snapshot(t)
            frust = (
                sg.get_frustration_index()
                if hasattr(sg, "get_frustration_index")
                else 0.0
            )
            frustrations.append((t, frust))
        return frustrations

    def evolve_signs_stochastic(
        self,
        t_start: float,
        t_end: float,
        dt: float,
        flip_rate: Optional[float] = None,
    ) -> None:
        """
        Evolve edge signs stochastically over a time interval.

        At each time step, edges flip with probability flip_rate * dt.

        Parameters
        ----------
        t_start : float
            Start time.
        t_end : float
            End time.
        dt : float
            Time step.
        flip_rate : float, optional
            Rate of sign flips per edge per unit time.
            Default: self.pflip.
        """
        if flip_rate is None:
            flip_rate = self.pflip

        t = t_start
        while t < t_end:
            G = self.get_snapshot(t)
            edges = list(G.edges())
            for u, v in edges:
                if self._temporal_rng.random() < flip_rate * dt:
                    old_sign = G[u][v].get("sign", 1)
                    self._sign_flips.append((t, u, v, -old_sign))
            self._times.add(t)
            t += dt

        self._invalidate_temporal_cache(t_start, t_end)

    # ------------------------------------------------------------------
    # Temporal metrics
    # ------------------------------------------------------------------

    def compute_temporal_metrics(self) -> dict:
        """
        Compute temporal network metrics.

        Returns
        -------
        dict
            Dictionary with temporal metrics:
            - "burstiness": Edge activity burstiness
            - "avg_duration": Average edge duration
            - "turnover": Edge turnover rate
        """
        durations = []
        for _, _, ts, te, _ in self._temporal_edges:
            if te != float("inf"):
                durations.append(te - ts)

        if not durations:
            return {"burstiness": 0.0, "avg_duration": float("inf"), "turnover": 0.0}

        durations = np.array(durations)
        avg = np.mean(durations)
        std = np.std(durations)

        burstiness = (std - avg) / (std + avg) if (std + avg) > 0 else 0.0

        t_range = self.time_range[1] - self.time_range[0]
        turnover = len(self._temporal_edges) / t_range if t_range > 0 else 0.0

        return {
            "burstiness": burstiness,
            "avg_duration": avg,
            "turnover": turnover,
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _merge_intervals(
        self, intervals: list[tuple[float, float]]
    ) -> list[tuple[float, float]]:
        """Merge overlapping time intervals."""
        if not intervals:
            return []
        sorted_intervals = sorted(intervals)
        merged = [sorted_intervals[0]]
        for start, end in sorted_intervals[1:]:
            if start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        return merged

    def _invalidate_temporal_cache(self, t_start: float, t_end: float) -> None:
        """Invalidate cached snapshots in a time range."""
        to_remove = [t for t in self._temporal_snapshots if t_start <= t < t_end]
        for t in to_remove:
            del self._temporal_snapshots[t]

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def from_edge_list(
        cls,
        edges: Sequence[tuple],
        n_nodes: Optional[int] = None,
        **kwargs,
    ) -> "TemporalGraphNX":
        """
        Create TemporalGraphNX from a list of temporal edges.

        Parameters
        ----------
        edges : sequence of tuples
            Each tuple is (u, v, t_start) or (u, v, t_start, t_end) or
            (u, v, t_start, t_end, attrs).
        n_nodes : int, optional
            Number of nodes.
        **kwargs
            Additional arguments for TemporalGraphNX.

        Returns
        -------
        TemporalGraphNX
            New temporal graph.
        """
        tg = cls(n_nodes=n_nodes, **kwargs)
        for edge in edges:
            u, v = edge[0], edge[1]
            t_start = edge[2] if len(edge) > 2 else 0.0
            t_end = edge[3] if len(edge) > 3 else None
            attrs = edge[4] if len(edge) > 4 else {}
            tg.add_temporal_edge(u, v, t_start=t_start, t_end=t_end, **attrs)
        return tg

    @classmethod
    def from_snapshot_sequence(
        cls,
        snapshots: Sequence[nx.Graph],
        times: Optional[Sequence[float]] = None,
        **kwargs,
    ) -> "TemporalGraphNX":
        """
        Create TemporalGraphNX from a sequence of snapshots.

        Parameters
        ----------
        snapshots : sequence of nx.Graph
            Graph snapshots in temporal order.
        times : sequence of float, optional
            Time points for each snapshot. Default: 0, 1, 2, ...
        **kwargs
            Additional arguments for TemporalGraphNX.

        Returns
        -------
        TemporalGraphNX
            New temporal graph.
        """
        if times is None:
            times = list(range(len(snapshots)))

        tg = cls(**kwargs)
        for t, G in zip(times, snapshots):
            tg.add_snapshot(t, graph=G)

        return tg
