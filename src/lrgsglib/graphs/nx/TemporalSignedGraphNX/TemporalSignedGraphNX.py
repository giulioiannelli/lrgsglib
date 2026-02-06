"""
TemporalSignedGraphNX: Temporal network with signed edges that can flip over time.

Extends TemporalGraphNX with specific support for signed edge dynamics,
including frustration evolution and sign flipping events.

Example
-------
>>> from lrgsglib.graphs.nx import TemporalSignedGraphNX
>>> tsg = TemporalSignedGraphNX(n_nodes=50, pflip=0.1, seed=42)
>>> tsg.add_edge(0, 1, t_start=0, sign=1)
>>> tsg.flip_random_edges(time=5.0, fraction=0.2)
>>> sg = tsg.get_signed_snapshot(time=6.0)
"""

from typing import Optional, Sequence
import numpy as np

from ..TemporalGraphNX.TemporalGraphNX import TemporalGraphNX
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


class TemporalSignedGraphNX(TemporalGraphNX):
    """
    Temporal network with signed edges that can flip over time.

    Extends TemporalGraphNX with specific support for signed edge dynamics,
    including frustration evolution and sign flipping events.

    Parameters
    ----------
    n_nodes : int, optional
        Number of nodes.
    pflip : float, default 0.0
        Probability of sign flip at each time step.
    seed : int, optional
        Random seed for reproducibility.
    **kwargs
        Additional arguments for TemporalGraphNX.

    Example
    -------
    >>> tsg = TemporalSignedGraphNX(n_nodes=50, pflip=0.1, seed=42)
    >>> tsg.add_edge(0, 1, t_start=0, sign=1)
    >>> tsg.flip_random_edges(time=5.0, fraction=0.2)
    >>> G_signed = tsg.get_signed_snapshot(time=6.0)
    """

    def __init__(
        self,
        n_nodes: Optional[int] = None,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ) -> None:
        super().__init__(n_nodes=n_nodes, **kwargs)
        self.pflip = pflip
        self.seed = seed
        self._rng = np.random.default_rng(seed)

        # Track sign flip events: list of (time, u, v, new_sign)
        self._sign_flips: list[tuple[float, int, int, int]] = []

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

        # Get current snapshot to know which edges exist
        G = self.get_snapshot(time)

        if edges is None:
            # Select random edges
            all_edges = list(G.edges())
            n_flip = int(len(all_edges) * fraction)
            flip_indices = self._rng.choice(len(all_edges), size=n_flip, replace=False)
            edges = [all_edges[i] for i in flip_indices]

        flipped = []
        for u, v in edges:
            if G.has_edge(u, v):
                old_sign = G[u][v].get("sign", 1)
                new_sign = -old_sign
                self._sign_flips.append((time, u, v, new_sign))
                flipped.append((u, v))

        self._times.add(time)
        self._invalidate_cache(time, float("inf"))

        return flipped

    def get_signed_snapshot(self, time: float) -> SignedGraphNX:
        """
        Get a SignedGraphNX snapshot at the given time.

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
                G[u][v]["weight"] = new_sign  # SignedGraphNX convention

        return SignedGraphNX(G)

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
            # Compute frustration index
            frust = sg.get_frustration_index() if hasattr(sg, "get_frustration_index") else 0.0
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
            # Get current edges
            G = self.get_snapshot(t)
            edges = list(G.edges())

            # Flip each edge with probability flip_rate * dt
            for u, v in edges:
                if self._rng.random() < flip_rate * dt:
                    old_sign = G[u][v].get("sign", 1)
                    self._sign_flips.append((t, u, v, -old_sign))

            self._times.add(t)
            t += dt

        self._invalidate_cache(t_start, t_end)
