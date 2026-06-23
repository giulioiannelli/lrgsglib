"""Negative-link (``nwDict``) pattern container for :class:`Lattice2DGT`.

Mirrors the key set of ``Lattice2DNX.nwContainer`` so notebooks can drive
identical frustration/defect placements on either engine
(``lat.nwDict[name]['G']`` → ``lat.flip_sel_edges(...)``):

``single`` / ``singleXERR`` / ``singleZERR`` — one defect at the central edge;
``rand`` / ``randXERR`` / ``randZERR`` — defects seeded on ``nflip`` random nodes.

The geometry itself (star, elementary cell, ball) is delegated to the
engine-neutral helpers in ``graphs/_shared/_nw_geometry.py``, which run over the
cross-engine ``get_graph_neighbors`` seam — so ZERR generalizes to whatever the
lattice's elementary face is (3-cycle triangular/kagome, 4-cycle square/octagon,
6-cycle honeycomb) without per-geometry code. Edges live in GT's integer node
space and are directly consumable by ``Lattice2DGT.flip_sel_edges``.
"""
from __future__ import annotations

import random
from typing import Any, List, Tuple

from ....config.const import COUNT_XERR_PATTERNS, L2D_ONREP
from ..._shared._nw_geometry import (
    elementary_cell_edges,
    star_edges,
)

Edge = Tuple[int, int]


class Lattice2DGTnwContainer(dict):
    """``nwDict`` for the graph-tool 2D lattice (see module docstring)."""

    def __init__(
        self,
        sg: "Any",
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
        self.centedge = {g: sg.get_central_edge(g) for g in self.rd}

        self["single"] = {g: [self.centedge[g]] for g in self.rd}
        self["singleZERR"] = {
            g: elementary_cell_edges(sg, self.centedge[g][0], g)
            for g in self.rd
        }
        self["singleXERR"] = {
            g: star_edges(sg, self.centedge[g][0], g) for g in self.rd
        }
        self["rand"] = {g: list(sg.fleset[g]) for g in self.rd}
        self["randZERR"] = {g: self._rand_pattern("ZERR", g) for g in self.rd}
        self["randXERR"] = {g: self._rand_pattern("XERR", g) for g in self.rd}

    def get_links_XERR(self, node: Any, on_g: str = L2D_ONREP) -> List[Edge]:
        """Star pattern (all edges incident to ``node``)."""
        return star_edges(self.sg, node, on_g)

    def get_links_ZERR(self, node: Any, on_g: str = L2D_ONREP) -> List[Edge]:
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


__all__ = ["Lattice2DGTnwContainer"]
