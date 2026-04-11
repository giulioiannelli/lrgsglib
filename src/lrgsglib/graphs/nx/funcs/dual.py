"""Planar dual graph construction for NetworkX graphs."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Optional, TYPE_CHECKING

import networkx as nx

if TYPE_CHECKING:
    from ..SignedGraphNX.SignedGraphNX import SignedGraphNX

__all__ = [
    "build_planar_dual",
    "build_signed_dual",
]


def build_planar_dual(
    G: nx.Graph,
    weight: str = "weight",
    include_outer_face: bool = True,
) -> tuple[nx.Graph, dict[str, Any]]:
    """Build the planar dual of a planar graph.

    Each face of the planar embedding becomes a node in the dual.
    Two dual nodes are connected if their faces share a primal edge,
    and the dual edge inherits the weight of the crossed primal edge.

    Parameters
    ----------
    G : nx.Graph
        A planar NetworkX graph with integer or hashable node labels.
    weight : str
        Edge attribute name to transfer to dual edges.
    include_outer_face : bool
        If True, include the infinite/outer face as a dual node.

    Returns
    -------
    dual : nx.Graph
        The dual graph. Nodes are integer face indices.
    face_info : dict
        ``'faces'``: list of face node-tuples,
        ``'edge_to_faces'``: dict mapping ``(u, v)`` to face-index pair,
        ``'outer_face_idx'``: index of the outer face (or None),
        ``'face_sizes'``: list of face perimeters.

    Raises
    ------
    nx.NetworkXException
        If *G* is not planar.
    """
    is_planar, embedding = nx.check_planarity(G)
    if not is_planar:
        raise nx.NetworkXException("Graph is not planar")

    # --- enumerate faces via half-edge traversal ---
    visited_half_edges: set[tuple] = set()
    faces: list[list] = []

    for v in embedding.nodes():
        for w in embedding.neighbors_cw_order(v):
            if (v, w) not in visited_half_edges:
                face = embedding.traverse_face(v, w, visited_half_edges)
                faces.append(face)

    # --- map each primal edge to its two bordering faces ---
    edge_to_faces: dict[tuple, list[int]] = defaultdict(list)

    for fi, face in enumerate(faces):
        n = len(face)
        for k in range(n):
            u, v = face[k], face[(k + 1) % n]
            e = (min(u, v), max(u, v))
            edge_to_faces[e].append(fi)

    # --- identify outer face (largest perimeter) ---
    outer_idx = max(range(len(faces)), key=lambda i: len(faces[i]))

    # --- build dual graph ---
    dual = nx.Graph()

    for fi in range(len(faces)):
        if not include_outer_face and fi == outer_idx:
            continue
        dual.add_node(fi, face_size=len(faces[fi]))

    for e, fp in edge_to_faces.items():
        if len(fp) != 2:
            continue
        f1, f2 = fp
        if not include_outer_face and (f1 == outer_idx or f2 == outer_idx):
            continue
        u, v = e
        w = G[u][v].get(weight, 1)
        dual.add_edge(f1, f2, weight=w)

    face_info = {
        "faces": [tuple(f) for f in faces],
        "edge_to_faces": dict(edge_to_faces),
        "outer_face_idx": outer_idx,
        "face_sizes": [len(f) for f in faces],
    }
    return dual, face_info


def build_signed_dual(
    sg: "SignedGraphNX",
    on_g: Optional[str] = None,
    include_outer_face: bool = True,
) -> "SignedGraphNX":
    """Build a :class:`SignedGraphNX` wrapping the planar dual.

    Edge signs transfer directly: a dual edge inherits the weight of
    the primal edge it crosses.  The returned object is ready for
    spectral analysis (``compute_k_eigvV``, ``compute_pinf``, etc.).

    Parameters
    ----------
    sg : SignedGraphNX
        Source signed graph (must be planar).
    on_g : str, optional
        Graph representation key. Defaults to ``sg.on_g``.
    include_outer_face : bool
        Whether to include the outer/infinite face.

    Returns
    -------
    SignedGraphNX
        A fresh SignedGraphNX with the dual topology and inherited signs.
    """
    from ..SignedGraphNX.SignedGraphNX import SignedGraphNX
    from ....config.const import SG_REPR

    if on_g is None:
        on_g = getattr(sg, "on_g", SG_REPR)

    G = sg.gr[on_g]
    dual_nx, face_info = build_planar_dual(G, include_outer_face=include_outer_face)

    # Relabel to consecutive integers for clean indexing
    dual_nx = nx.convert_node_labels_to_integers(dual_nx)

    # Wrap in SignedGraphNX with pflip=0 (signs already baked in)
    dual_sg = SignedGraphNX(
        dual_nx,
        pflip=0.0,
        make_dir_tree=False,
    )

    # Fix edge sets to reflect actual negative weights from inheritance
    dual_sg.upd_edge_sets(on_g=dual_sg.on_g)
    dual_sg.pflip = (
        len(dual_sg.fleset[dual_sg.on_g]) / dual_sg.Ne
        if dual_sg.Ne > 0
        else 0.0
    )
    dual_sg.Ne_flips = len(dual_sg.fleset[dual_sg.on_g])

    # Store provenance metadata
    dual_sg._dual_face_info = face_info
    dual_sg._dual_of = sg

    return dual_sg
