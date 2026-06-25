"""Engine-neutral negative-link (``nwDict``) pattern container.

A single base, ``NwContainer``, builds the negative-link defect patterns for
**any** signed graph on **either** engine (NX or GT). It is the default
``nwContainer`` on both ``SignedGraph`` bases, so every graph — lattices, ER, BA,
Holme-Kim, … — gets the geometry-free patterns (``rand``, ``randXERR``) for free.
Two class-level flags add the geometry-dependent patterns where they apply:

==================  ==============  =============  =================================
graph class         build_single    build_zerr     keys built
==================  ==============  =============  =================================
default (ER/BA/…)   False           False          rand, randXERR
Lattice3D           True            False          single, singleXERR, rand, randXERR
Lattice2D           True            True           + singleZERR, randZERR
==================  ==============  =============  =================================

The geometry is delegated to the engine-neutral helpers in ``_nw_geometry``
(star = XERR, elementary cell = generalized ZERR, all over the cross-engine
``get_graph_neighbors`` seam), so there is no per-engine or per-geometry pattern
code here. ``single*`` patterns require the graph to implement
``get_central_edge`` (lattices do; non-geometric graphs do not, hence
``build_single`` defaults ``False``). ``GTnwContainer`` is kept as a back-compat
alias.
"""
from __future__ import annotations

import random
from typing import Any, List, Tuple

import numpy as np

from ...config.const import COUNT_XERR_PATTERNS, SG_REPR
from ._nw_geometry import star_edges

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


class NwContainer(dict):
    """Engine-neutral negative-link pattern dict (see module doc)."""

    #: build ``single`` / ``singleXERR`` (needs ``sg.get_central_edge``).
    build_single: bool = False
    #: build the elementary-cell ZERR patterns (``singleZERR`` / ``randZERR``).
    build_zerr: bool = False

    def __init__(
        self,
        sg: Any,
        iterable=(),
        constant: Any = None,
        *,
        build_single: bool | None = None,
        build_zerr: bool | None = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.update((key, constant) for key in iterable)

        # Per-build overrides win over the class-level defaults. Lattices set
        # the flags True as class attrs; geometry-free graphs (ER/BA/…) leave
        # them False but the engine base passes an override here when a
        # structured ``single*`` / ``*ZERR`` disorder support is requested, so
        # those patterns build on demand for *any* subclass.
        bs = self.build_single if build_single is None else build_single
        bz = self.build_zerr if build_zerr is None else build_zerr

        self.sg = sg
        self.rd = sg.graph_reprs
        # Sample the random seed nodes per representation: a graph may label its
        # reprs differently (e.g. a lattice's ``G`` are ints, ``H`` are coordinate
        # tuples), so the seeds for the ``H`` patterns must come from ``H``'s own
        # node set, not the default repr's.
        self.rNodeFlip = {
            g: random.sample(sg.get_nodes_list(g), sg.nflip) for g in self.rd
        }

        if bs:
            self.centedge = {g: sg.get_central_edge(g) for g in self.rd}
            self["single"] = {g: [self.centedge[g]] for g in self.rd}
            self["singleXERR"] = {
                g: star_edges(sg, self.centedge[g][0], g) for g in self.rd
            }
            if bz:
                self["singleZERR"] = {
                    g: sg.cell_edges(self.centedge[g][0], g)
                    for g in self.rd
                }

        self["rand"] = {g: list(sg.fleset[g]) for g in self.rd}
        self["randXERR"] = {g: self._rand_pattern("XERR", g) for g in self.rd}
        if bz:
            self["randZERR"] = {g: self._rand_pattern("ZERR", g) for g in self.rd}

    def get_links_XERR(self, node: Any, on_g: str = SG_REPR) -> List[Edge]:
        """Star pattern (all edges incident to ``node``)."""
        return star_edges(self.sg, node, on_g)

    def get_links_ZERR(self, node: Any, on_g: str = SG_REPR) -> List[Edge]:
        """Elementary-cell pattern (smallest cycle through ``node``)."""
        return self.sg.cell_edges(node, on_g)

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
                    for e in self.sg.cell_edges(i, on_g)
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


# Back-compat alias (was GT-only before the container became engine-neutral).
GTnwContainer = NwContainer

__all__ = ["NwContainer", "GTnwContainer", "geometric_central_edge"]
