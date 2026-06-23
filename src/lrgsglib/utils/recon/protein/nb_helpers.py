"""High-level helpers consumed by the CHL-SF04 notebook.

Hosted here (rather than as inline ``def _foo`` cells) so the report file
stays compute-light per ``.agents/rules/notebook-hygiene.md``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import networkx as nx
import numpy as np

from ....config.const import CHL_SF04_N_FEATURES
from .feature_extraction import assign_secondary_structure_from_coords
from .metrics import kabsch_align, kabsch_rmsd, q3_score
from .substrate import SubstrateSpec, build_substrate_basis

__all__ = [
    "load_protein_corpus",
    "load_tmd_bundles",
    "interp_mse_curve",
    "mse_per_substrate",
    "pflip_curves",
    "auc_per_substrate",
    "q3_emergence_per_substrate",
    "q3_curve_per_substrate",
    "ss_segments",
    "layout_positions_2d",
    "render_substrate_panel",
    "reconstruct_protein",
]


# ---------------------------------------------------------------------------
# Bundle / corpus loaders
# ---------------------------------------------------------------------------


def load_protein_corpus(
    corpus_file: Path,
    pdb_dir: Path,
    download_pdb_fn,
    extract_ca_fn,
) -> List[Dict[str, Any]]:
    """Load the frozen corpus index and re-extract centred CA coords + SS."""
    corpus = json.loads(corpus_file.read_text())
    proteins: List[Dict[str, Any]] = []
    for entry in corpus:
        pdb_id = entry["pdb_id"]
        coords = extract_ca_fn(download_pdb_fn(pdb_id, pdb_dir))
        ss = assign_secondary_structure_from_coords(coords)
        proteins.append({
            "pdb_id": pdb_id,
            "n_residues": int(len(coords)),
            "coords": coords - coords.mean(axis=0),
            "ss": ss,
        })
    return proteins


def load_tmd_bundles(tmd_dir: Path) -> List[Dict[str, Any]]:
    """Glob ``*.npz`` bundles + their sibling ``*.json`` attribute files."""
    bundles: List[Dict[str, Any]] = []
    for npz_path in sorted(tmd_dir.glob("*.npz")):
        attrs = json.loads(npz_path.with_suffix(".json").read_text())
        bundles.append({
            "npz": npz_path,
            "attrs": attrs,
            "label": attrs["label"],
            "spec": attrs["spec"],
            "pdb_id": attrs["pdb_id"],
            "N": attrs["N"],
        })
    return bundles


# ---------------------------------------------------------------------------
# MSE-curve helpers
# ---------------------------------------------------------------------------


def interp_mse_curve(
    npz_path: Path,
    k_grid: np.ndarray,
    floor: float = 0.0,
) -> np.ndarray:
    """Re-sample a stored ``mse_curve`` onto a common ``k/N`` grid, normalised by MSE(0).

    The stored curve is indexed by ``k = 0 .. N-1`` (the trivial ``k = N``
    full-reconstruction point is *not* persisted; it would only force every
    log-y curve to zero). ``floor`` clips the normalised result from below so
    the log-scale plot stays readable.
    """
    arr = np.load(npz_path)
    mse = np.asarray(arr["mse_curve"], dtype=np.float64)
    k_native = np.arange(len(mse)) / max(len(mse), 1)  # k/N for k = 0..N-1
    norm = mse / max(mse[0], 1e-12)
    if floor > 0:
        norm = np.maximum(norm, floor)
    return np.interp(k_grid, k_native, norm)


def mse_per_substrate(
    bundles: Sequence[Dict[str, Any]],
    k_grid: np.ndarray,
    pflip: float = 0.0,
    floor: float = 0.0,
) -> Dict[str, np.ndarray]:
    """Group MSE curves by substrate label, optionally filtering by ``pflip``."""
    out: Dict[str, List[np.ndarray]] = {}
    for b in bundles:
        if float(b["spec"]["pflip"]) != pflip:
            continue
        out.setdefault(b["label"], []).append(
            interp_mse_curve(b["npz"], k_grid, floor=floor)
        )
    return {label: np.array(curves) for label, curves in out.items() if curves}


def _matches_pflip_family(label: str, family_label: str) -> bool:
    """A bundle ``label`` matches ``family_label`` (which may carry a ``_pfX`` tag)."""
    return label.split("_pf")[0] == family_label.split("_pf")[0]


def pflip_curves(
    bundles: Sequence[Dict[str, Any]],
    family_label: str,
    pflip_grid: Sequence[float],
    k_grid: np.ndarray,
    floor: float = 0.0,
) -> Dict[float, np.ndarray]:
    """For one substrate family, gather MSE curves at each ``pflip``."""
    out: Dict[float, List[np.ndarray]] = {p: [] for p in pflip_grid}
    for b in bundles:
        if not _matches_pflip_family(b["label"], family_label):
            continue
        pf = float(b["spec"]["pflip"])
        if pf not in out:
            continue
        out[pf].append(interp_mse_curve(b["npz"], k_grid, floor=floor))
    return {p: np.array(v) for p, v in out.items() if v}


def auc_per_substrate(
    mse_groups: Dict[str, np.ndarray],
    k_grid: np.ndarray,
    k_lo: float = 0.05,
    k_hi: float = 0.50,
) -> Dict[str, float]:
    """Trapezoidal area-under-mean-MSE in ``[k_lo, k_hi]``; smaller = better."""
    mask = (k_grid >= k_lo) & (k_grid <= k_hi)
    return {
        label: float(np.trapz(curves.mean(axis=0)[mask], k_grid[mask]))
        for label, curves in mse_groups.items()
    }


# ---------------------------------------------------------------------------
# Q3 / SS helpers
# ---------------------------------------------------------------------------


def reconstruct_protein(
    protein: Dict[str, Any],
    spec: SubstrateSpec,
    k_frac: float,
    spectrum_backend: str = "cupy",
    n_features: int = CHL_SF04_N_FEATURES,
) -> Tuple[np.ndarray, float]:
    """Rebuild basis on-the-fly, project, truncate to ``k_frac`` modes, return (coords, RMSD).

    The substrate is sized so ``N >= n_features`` (default 1024) regardless of
    protein length, matching the lab compute pass; the protein coordinates
    occupy the first ``3*n_res`` entries and are zero-padded to ``N``.
    """
    n_res = protein["n_residues"]
    n_dof = 3 * n_res
    basis, sg, _ = build_substrate_basis(
        spec, n_target=n_features, spectrum_backend=spectrum_backend
    )
    if n_dof > sg.N:
        raise RuntimeError(
            f"{protein['pdb_id']}: 3*n_residues={n_dof} > N={sg.N}; raise n_features."
        )
    x_padded = np.zeros(sg.N, dtype=np.float64)
    x_padded[:n_dof] = protein["coords"].flatten()
    k = max(1, int(round(k_frac * sg.N)))
    recon_full = basis[:k].T @ (basis[:k] @ x_padded)
    recon = recon_full[:n_dof].reshape(n_res, 3)
    return kabsch_align(recon, protein["coords"]), kabsch_rmsd(recon, protein["coords"])


def q3_curve_per_substrate(
    proteins: Sequence[Dict[str, Any]],
    bundles: Sequence[Dict[str, Any]],
    label: str,
    k_frac_grid: Sequence[float],
    spectrum_backend: str = "cupy",
    n_features: int = CHL_SF04_N_FEATURES,
) -> np.ndarray:
    """Q3 vs ``k/N`` for one substrate, averaged later by the caller."""
    matrix: List[List[float]] = []
    for protein in proteins:
        spec_dict = next(
            (b["spec"] for b in bundles
             if b["pdb_id"] == protein["pdb_id"] and b["label"] == label
             and float(b["spec"]["pflip"]) == 0.0),
            None,
        )
        if spec_dict is None:
            continue
        spec = SubstrateSpec(**{**spec_dict, "extra_kw": tuple()})
        n_res = protein["n_residues"]
        n_dof = 3 * n_res
        basis, sg, _ = build_substrate_basis(
            spec, n_target=n_features, spectrum_backend=spectrum_backend
        )
        x_padded = np.zeros(sg.N, dtype=np.float64)
        x_padded[:n_dof] = protein["coords"].flatten()
        row: List[float] = []
        for k_frac in k_frac_grid:
            k = max(1, int(round(k_frac * sg.N)))
            recon = (basis[:k].T @ (basis[:k] @ x_padded))[:n_dof].reshape(n_res, 3)
            ss_pred = assign_secondary_structure_from_coords(recon)
            row.append(q3_score(ss_pred, protein["ss"]))
        matrix.append(row)
    return np.array(matrix) if matrix else np.empty((0, len(k_frac_grid)))


def q3_emergence_per_substrate(
    proteins: Sequence[Dict[str, Any]],
    bundles: Sequence[Dict[str, Any]],
    label: str,
    k_frac_grid: Sequence[float],
    threshold: float,
    spectrum_backend: str = "cupy",
) -> float:
    """Mean (over corpus) of the smallest ``k/N`` whose Q3 ≥ ``threshold``."""
    curve = q3_curve_per_substrate(
        proteins, bundles, label, k_frac_grid, spectrum_backend=spectrum_backend
    )
    if curve.size == 0:
        return 1.0
    k_arr = np.asarray(k_frac_grid, dtype=np.float64)
    emergence: List[float] = []
    for row in curve:
        hits = np.where(row >= threshold)[0]
        emergence.append(float(k_arr[hits[0]]) if hits.size else 1.0)
    return float(np.mean(emergence)) if emergence else 1.0


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def ss_segments(coords: np.ndarray, ss: Sequence[str], ss_colours: Dict[str, str]):
    """Backbone segments + per-segment colours for a Line3DCollection."""
    pts = np.asarray(coords)
    segs = [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    cols = [ss_colours.get(ss[i], "#888888") for i in range(len(pts) - 1)]
    return segs, cols


def layout_positions_2d(sg, fallback_seed: int = 0) -> Dict[Any, np.ndarray]:
    """Return a 2D layout dict for a SignedGraph, preferring stored ``pos`` attrs."""
    H = sg.gr.get("H") if isinstance(sg.gr, dict) else None
    if H is not None:
        pos = nx.get_node_attributes(H, "pos")
        if pos:
            return {n: np.asarray(p)[:2] for n, p in pos.items()}
    G = sg.gr["G"] if isinstance(sg.gr, dict) else getattr(sg, "G")
    return nx.spring_layout(G, seed=fallback_seed)


def render_substrate_panel(
    ax,
    spec: SubstrateSpec,
    n_target: int,
    title: str,
    spectrum_backend: str = "cupy",
) -> None:
    """Render one panel of the substrate-zoo figure (Fig. 1).

    Edge colour: red for negative, grey for positive. Imports
    ``LineCollection`` lazily so the helper module stays import-light.
    """
    from matplotlib.collections import LineCollection

    _, sg, _ = build_substrate_basis(
        spec, n_target=n_target, spectrum_backend=spectrum_backend
    )
    G = sg.gr["G"] if isinstance(sg.gr, dict) else getattr(sg, "G")
    pos = layout_positions_2d(sg)
    edge_pos: List[Tuple[np.ndarray, np.ndarray]] = []
    edge_neg: List[Tuple[np.ndarray, np.ndarray]] = []
    for u, v, d in G.edges(data=True):
        seg = (pos[u], pos[v])
        (edge_neg if d.get("weight", 1) < 0 else edge_pos).append(seg)
    if edge_pos:
        ax.add_collection(LineCollection(
            edge_pos, colors="0.6", linewidths=0.4, alpha=0.7))
    if edge_neg:
        ax.add_collection(LineCollection(
            edge_neg, colors="#d62728", linewidths=0.6, alpha=0.9))
    pts = np.array([pos[n] for n in G.nodes])
    ax.scatter(pts[:, 0], pts[:, 1], s=4, c="#222222", zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=9)
    ax.autoscale()
