"""Engine-agnostic 2D lattice drawing.

A single :func:`draw` routine bound as a method on both ``Lattice2DNX`` and
``Lattice2DGT``.  It relies only on accessors present on *both* engines
(:meth:`get_node_attributes`, :meth:`get_edges_with_weights`), so identical
user code works whatever the active backend::

    lat = Lattice2D(side1=24, geo='sqr', pbc=True, with_positions=True,
                    engine='nx')   # or engine='gt' -- same call below
    lat.draw()

Rendering goes straight through matplotlib (a single ``LineCollection`` for the
bonds plus one scatter for the sites) instead of graph-tool's Cairo canvas,
which is fast, gives a true 2D layout, and behaves the same on every engine.
"""
from __future__ import annotations

from typing import Any, Optional, Tuple

import numpy as np


def draw(
    self,
    ax: Optional[Any] = None,
    *,
    node_size: float = 20.0,
    node_color: Any = "#333333",
    edge_color: Optional[Any] = None,
    width: float = 0.8,
    pec: Any = "blue",
    nec: Any = "red",
    hide_wrap: bool = False,
    wrap_factor: float = 1.8,
    with_labels: bool = False,
    figsize: Tuple[float, float] = (6.0, 6.0),
    **kwargs: Any,
):
    """Draw the 2D lattice using its stored node positions.

    Identical behaviour on the NetworkX and graph-tool backends.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes. A new ``(figure, axes)`` pair is created when omitted.
    node_size : float
        Marker area for sites (matplotlib ``scatter`` ``s``).
    node_color : color
        Site marker colour.
    edge_color : color, optional
        If given, every bond uses this single colour (mirrors ``nx.draw``).
        If ``None`` (default), bonds are coloured by sign with *pec*/*nec*.
    width : float
        Bond line width.
    pec, nec : color
        Positive / negative bond colours when colouring by sign.
    hide_wrap : bool
        If True, omit periodic wrap-around bonds (those much longer than a
        typical bond), giving a clean open-lattice look on a torus.
    wrap_factor : float
        A bond counts as a wrap bond when its length exceeds ``wrap_factor``
        times the median bond length.
    with_labels : bool
        Annotate each site with its integer index.
    figsize : (float, float)
        Figure size used only when *ax* is None.
    **kwargs
        Forwarded to ``matplotlib.collections.LineCollection``.

    Returns
    -------
    matplotlib.axes.Axes
        The axes the lattice was drawn on.

    Raises
    ------
    ValueError
        If no node positions are stored. Construct the lattice with
        ``with_positions=True``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    pos = self.get_node_attributes("pos")
    if not pos:
        raise ValueError(
            "No node positions available to draw. Construct the lattice "
            "with `with_positions=True`."
        )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    # Build bond segments + per-bond sign using engine-agnostic accessors.
    segments, lengths, signs = [], [], []
    for u, v, w in self.get_edges_with_weights():
        if u not in pos or v not in pos:
            continue
        p, q = np.asarray(pos[u], dtype=float), np.asarray(pos[v], dtype=float)
        segments.append((p, q))
        lengths.append(float(np.hypot(*(q - p))))
        signs.append(w)

    lengths = np.asarray(lengths)
    keep = np.ones(len(segments), dtype=bool)
    if hide_wrap and lengths.size:
        med = np.median(lengths)
        if med > 0:
            keep = lengths <= wrap_factor * med

    kept_segments = [s for s, k in zip(segments, keep) if k]
    if edge_color is not None:
        colors: Any = edge_color
    else:
        colors = [nec if s < 0 else pec for s, k in zip(signs, keep) if k]

    ax.add_collection(
        LineCollection(kept_segments, colors=colors, linewidths=width, **kwargs)
    )

    coords = np.asarray([pos[n] for n in pos], dtype=float)
    ax.scatter(coords[:, 0], coords[:, 1], s=node_size, c=node_color, zorder=2)

    if with_labels:
        for n, xy in pos.items():
            ax.annotate(str(n), (float(xy[0]), float(xy[1])),
                        fontsize=6, ha="center", va="center")

    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.margins(0.05)
    ax.set_axis_off()
    return ax
