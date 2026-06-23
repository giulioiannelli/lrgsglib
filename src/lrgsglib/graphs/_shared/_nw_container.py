"""Configurable negative-link (``nwDict``) pattern container for GT graphs.

A single base, ``GTnwContainer``, builds the negative-link defect patterns that
the NX engine exposes per graph class, driven by two class-level flags so each
GT graph type opts into exactly the patterns its NX counterpart has:

==================  ==============  =============  =================================
graph class         build_single    build_zerr     keys built
==================  ==============  =============  =================================
ErdosRenyiGT        False           False          rand, randXERR
Lattice3DGT         True            False          single, singleXERR, rand, randXERR
Lattice2DGT         True            True           + singleZERR, randZERR
==================  ==============  =============  =================================

The geometry is delegated to the engine-neutral helpers in ``_nw_geometry``
(star = XERR, elementary cell = generalized ZERR, all over the cross-engine
``get_graph_neighbors`` seam), so there is no per-engine or per-geometry pattern
code here. ``single*`` patterns require the graph to implement
``get_central_edge`` (lattices do; ErdosRenyi does not, hence ``build_single``).
"""
from __future__ import annotations

import random
from typing import Any, List, Tuple

import numpy as np

from ...config.const import COUNT_XERR_PATTERNS, SG_REPR
from ._nw_geometry import elementary_cell_edges, star_edges

Edge = Tuple[int, int]


def geometric_central_edge(sg: Any, on_g: str = SG_REPR) -> Edge:
    """A bulk edge nearest the geometric centre of a GT lattice (any dimension).

    Operates in GT's integer node space using the stored ``_pos`` vertex
    property: the chosen edge is the incident edge of the centroid-nearest
    vertex whose midpoint lies closest to the centroid. Deterministic and
    translation-stable. Shared by ``Lattice2DGT`` and ``Lattice3DGT`` (their
    ``pos`` arrays are 2D / 3D respectively; the logic is dimension-agnostic).
    ``on_g`` is accepted for API parity (GT has a single representation).
    """
    pos = np.array([list(sg._pos[sg.G.vertex(i)]) for i in range(sg.N)])
    centroid = pos.mean(axis=0)
    central = int(np.argmin(((pos - centroid) ** 2).sum(axis=1)))
    nbrs = sg.get_graph_neighbors(central)
    if not nbrs:
        raise ValueError(
            "central node has no neighbours; cannot pick a central edge"
        )
    nb = min(
        nbrs,
        key=lambda j: ((0.5 * (pos[central] + pos[j]) - centroid) ** 2).sum(),
    )
    return (central, nb)


class GTnwContainer(dict):
    """Negative-link pattern dict for the graph-tool engine (see module doc)."""

    #: build ``single`` / ``singleXERR`` (needs ``sg.get_central_edge``).
    build_single: bool = False
    #: build the elementary-cell ZERR patterns (``singleZERR`` / ``randZERR``).
    build_zerr: bool = False

    def __init__(
        self,
        sg: Any,
        iterable=(),
        constant: Any = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.update((key, constant) for key in iterable)

        self.sg = sg
        self.rd = sg.graph_reprs
        self.rNodeFlip = {
            g: random.sample(sg.get_nodes_list(), sg.nflip) for g in self.rd
        }

        if self.build_single:
            self.centedge = {g: sg.get_central_edge(g) for g in self.rd}
            self["single"] = {g: [self.centedge[g]] for g in self.rd}
            self["singleXERR"] = {
                g: star_edges(sg, self.centedge[g][0], g) for g in self.rd
            }
            if self.build_zerr:
                self["singleZERR"] = {
                    g: elementary_cell_edges(sg, self.centedge[g][0], g)
                    for g in self.rd
                }

        self["rand"] = {g: list(sg.fleset[g]) for g in self.rd}
        self["randXERR"] = {g: self._rand_pattern("XERR", g) for g in self.rd}
        if self.build_zerr:
            self["randZERR"] = {g: self._rand_pattern("ZERR", g) for g in self.rd}

    def get_links_XERR(self, node: Any, on_g: str = SG_REPR) -> List[Edge]:
        """Star pattern (all edges incident to ``node``)."""
        return star_edges(self.sg, node, on_g)

    def get_links_ZERR(self, node: Any, on_g: str = SG_REPR) -> List[Edge]:
        """Elementary-cell pattern (smallest cycle through ``node``)."""
        return elementary_cell_edges(self.sg, node, on_g)

    def _rand_pattern(self, mode: str, on_g: str) -> List[Edge]:
        nodes = self.rNodeFlip[on_g]
        match mode:
            case "XERR":
                if COUNT_XERR_PATTERNS:
                    pattern = [
                        e for i in nodes for e in star_edges(self.sg, i, on_g)
                    ]
                else:
                    pattern = self._filtered_xerr(on_g)
            case "ZERR":
                pattern = [
                    e
                    for i in nodes
                    for e in elementary_cell_edges(self.sg, i, on_g)
                ]
            case _:
                pattern = []
        return list(set(pattern))

    def _filtered_xerr(self, on_g: str) -> List[Edge]:
        """XERR with the NX filter: skip seed nodes that have a neighbour whose
        every incident edge is already negative (mirrors the NX container)."""
        tmp = list(self.rNodeFlip[on_g])
        idx = 0
        pattern: List[Edge] = []
        while idx < len(tmp):
            node = tmp[idx]
            has_all_neg_neighbor = any(
                all(
                    w == -1
                    for _, w in self.sg.get_neighbors_with_weights(nn)
                )
                for nn in self.sg.get_graph_neighbors(node, on_g)
            )
            if has_all_neg_neighbor:
                tmp.pop(idx)
            else:
                pattern.extend(star_edges(self.sg, node, on_g))
                idx += 1
        return pattern


__all__ = ["GTnwContainer"]
