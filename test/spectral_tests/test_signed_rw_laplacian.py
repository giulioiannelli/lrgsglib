"""
Tests for the optional signed random-walk / symmetric-normalized Laplacians.

Covers the three ``laplacian_type`` values threaded through the SignedGraph
spectral machinery (Kunegis 2010, kunegis2010spectral.pdf §3.2-3.4):

    'signed' -> L     = D_s - A                  (combinatorial, symmetric, default)
    'sym'    -> L_sym = I - D_s^-1/2 A D_s^-1/2  (symmetric normalized)
    'rw'     -> L_rw  = I - D_s^-1 A             (random walk, NON-symmetric)

The default ('signed') path must be unchanged. ``rw`` and ``sym`` are distinct
matrices/solvers/eigenvectors but isospectral (real eigenvalues in [0, 2]).
"""

import numpy as np
import networkx as nx
import pytest

from lrgsglib.graphs import SignedGraph
from lrgsglib.graphs.nx.funcs.spectral import (
    signed_laplacian_matrix,
    signed_sym_laplacian_matrix,
    signed_rw_laplacian_matrix,
)
from lrgsglib.utils.lrg.spectral import get_graph_lspectrum, get_graph_lspectrum_rw
from lrgsglib.config.const import (
    SG_LAPL_SIGNED,
    SG_LAPL_RW,
    SG_LAPL_SYM,
    SG_LAPL_DEFAULT_TYPE,
    SG_LAPL_RW_IMAG_TOL,
)


# --------------------------------------------------------------------------- #
# Small hand-verifiable graphs
# --------------------------------------------------------------------------- #
def _weighted(g):
    nx.set_edge_attributes(g, 1.0, "weight")
    return g


def _signed_graph(g):
    """Wrap a weighted nx.Graph as a NetworkX SignedGraph."""
    return SignedGraph(G=g, engine="nx")


def _balanced_path(n=7):
    """Path graph (a tree -> always balanced). Degrees 1,2,...,2,1 (irregular)."""
    return _signed_graph(_weighted(nx.path_graph(n)))


def _frustrated_triangle():
    """3-cycle with exactly one negative edge -> unbalanced (frustrated)."""
    g = _weighted(nx.cycle_graph(3))
    g[0][1]["weight"] = -1.0
    return _signed_graph(g)


def _irregular_signed():
    """Star + tail: degrees vary strongly, so L_rw != L_sym (non-symmetric)."""
    g = nx.star_graph(6)       # center 0 (deg 6), leaves deg 1
    g.add_edge(6, 7)
    g.add_edge(7, 8)
    g = _weighted(g)
    g[0][1]["weight"] = -1.0   # a sign, still a tree -> balanced
    return _signed_graph(g)


# --------------------------------------------------------------------------- #
# 1. Default 'signed' path is unchanged
# --------------------------------------------------------------------------- #
def test_default_signed_path_unchanged():
    sg = _irregular_signed()
    sg.compute_laplacian_spectrum()                     # default
    ev_default = np.sort(sg.eigv.copy())

    sg2 = _irregular_signed()
    sg2.compute_laplacian_spectrum(laplacian_type=SG_LAPL_SIGNED)
    np.testing.assert_allclose(np.sort(sg2.eigv), ev_default, atol=1e-10)

    # ... and equals the combinatorial signed Laplacian free function
    L = signed_laplacian_matrix(sg.gr[sg.on_g]).toarray()
    np.testing.assert_allclose(
        np.sort(np.linalg.eigvalsh(L)), ev_default, atol=1e-9
    )


# --------------------------------------------------------------------------- #
# 2. Non-degeneracy: rw is genuinely non-symmetric, sym is symmetric
# --------------------------------------------------------------------------- #
def test_rw_nonsymmetric_sym_symmetric():
    sg = _irregular_signed()
    Lrw = sg.get_signed_rw_laplacian(sym=False).toarray()
    Lsym = sg.get_signed_rw_laplacian(sym=True).toarray()
    assert not np.allclose(Lrw, Lrw.T), "L_rw must be non-symmetric on an irregular graph"
    assert np.allclose(Lsym, Lsym.T), "L_sym must be symmetric"
    # rw and sym are different matrices
    assert not np.allclose(Lrw, Lsym)


# --------------------------------------------------------------------------- #
# 3. Isospectrality: eigvals(rw, general) == eigvals(sym, symmetric)
# --------------------------------------------------------------------------- #
def test_isospectral_rw_equals_sym():
    sg = _frustrated_triangle()
    Lrw = sg.get_signed_rw_laplacian(sym=False).toarray()
    Lsym = sg.get_signed_rw_laplacian(sym=True).toarray()
    ev_rw = np.sort(np.linalg.eigvals(Lrw).real)
    ev_sym = np.sort(np.linalg.eigvalsh(Lsym))
    np.testing.assert_allclose(ev_rw, ev_sym, atol=1e-9)


# --------------------------------------------------------------------------- #
# 4. rw eigenvalues are real (through the compute method's real-cast path)
# --------------------------------------------------------------------------- #
def test_rw_eigenvalues_real_and_bounded():
    sg = _irregular_signed()
    sg.compute_laplacian_spectrum_weigV(laplacian_type=SG_LAPL_RW)
    assert not np.iscomplexobj(sg.eigv)
    assert np.all(np.isfinite(sg.eigv))
    assert sg.eigv.min() >= -1e-8 and sg.eigv.max() <= 2.0 + 1e-8


# --------------------------------------------------------------------------- #
# 5. sym eigenvectors are orthonormal; rw eigenvectors are not
# --------------------------------------------------------------------------- #
def test_eigenvectors_sym_orthonormal_rw_not():
    sg = _irregular_signed()
    sg.compute_laplacian_spectrum_weigV(
        laplacian_type=SG_LAPL_SYM, transpose=False, flip_to_pos=False
    )
    Vsym = sg.eigV
    assert np.allclose(Vsym.T @ Vsym, np.eye(Vsym.shape[1]), atol=1e-8)

    sg2 = _irregular_signed()
    sg2.compute_laplacian_spectrum_weigV(
        laplacian_type=SG_LAPL_RW, transpose=False, flip_to_pos=False
    )
    Vrw = sg2.eigV
    assert not np.allclose(Vrw.T @ Vrw, np.eye(Vrw.shape[1]), atol=1e-6)


# --------------------------------------------------------------------------- #
# 6. Builder method matches the canonical free-function kernels
# --------------------------------------------------------------------------- #
def test_method_matches_free_function():
    sg = _irregular_signed()
    G = sg.gr[sg.on_g]
    np.testing.assert_allclose(
        sg.get_signed_rw_laplacian(sym=True).toarray(),
        signed_sym_laplacian_matrix(G).toarray(),
        atol=1e-12,
    )
    np.testing.assert_allclose(
        sg.get_signed_rw_laplacian(sym=False).toarray(),
        signed_rw_laplacian_matrix(G).toarray(),
        atol=1e-12,
    )


# --------------------------------------------------------------------------- #
# 7. Balanced graph -> smallest normalized eigenvalue is 0 (Kunegis Thm 4.1)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("lap_type", [SG_LAPL_SYM, SG_LAPL_RW])
def test_balanced_graph_zero_eigenvalue(lap_type):
    sg = _balanced_path(8)  # a tree is always balanced
    sg.compute_laplacian_spectrum(laplacian_type=lap_type)
    ev = np.sort(sg.eigv)
    assert ev[0] == pytest.approx(0.0, abs=1e-9)
    assert ev[-1] <= 2.0 + 1e-9


# --------------------------------------------------------------------------- #
# 8. Unbalanced graph -> smallest normalized eigenvalue is strictly > 0 (Thm 4.4)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("lap_type", [SG_LAPL_SYM, SG_LAPL_RW])
def test_unbalanced_graph_positive_eigenvalue(lap_type):
    sg = _frustrated_triangle()
    sg.compute_laplacian_spectrum(laplacian_type=lap_type)
    assert np.sort(sg.eigv)[0] > 1e-9


# --------------------------------------------------------------------------- #
# 9. Spectrum cache is laplacian_type-aware
# --------------------------------------------------------------------------- #
def test_spectrum_cache_is_type_aware():
    sg = _irregular_signed()
    sg.compute_laplacian_spectrum(laplacian_type=SG_LAPL_SIGNED)
    ev_signed = np.sort(sg.eigv.copy())
    assert sg._spectrum_laplacian_type == SG_LAPL_SIGNED

    sg.compute_laplacian_spectrum(laplacian_type=SG_LAPL_SYM)
    ev_sym = np.sort(sg.eigv.copy())
    assert sg._spectrum_laplacian_type == SG_LAPL_SYM
    # combinatorial vs normalized differ (different operators)
    assert not np.allclose(ev_signed, ev_sym)

    # re-request signed -> cache must return the signed spectrum, not the sym one
    sg.compute_laplacian_spectrum(laplacian_type=SG_LAPL_SIGNED)
    np.testing.assert_allclose(np.sort(sg.eigv), ev_signed, atol=1e-10)


# --------------------------------------------------------------------------- #
# 10. Zero-degree (isolated) node -> finite matrices, eigenvalue 1, no NaN
# --------------------------------------------------------------------------- #
def test_isolated_node_no_nan():
    g = nx.Graph()  # noqa: raw-graph
    g.add_nodes_from([0, 1, 2])
    g.add_edge(0, 1, weight=1.0)        # node 2 is isolated
    sg = _signed_graph(g)
    for sym in (True, False):
        L = sg.get_signed_rw_laplacian(sym=sym).toarray()
        assert np.all(np.isfinite(L))
        ev = np.sort(np.linalg.eigvals(L).real)
        assert np.any(np.isclose(ev, 1.0)), "isolated node should contribute eigenvalue 1"


# --------------------------------------------------------------------------- #
# 11. Spectral gap order parameter is type-aware
# --------------------------------------------------------------------------- #
def test_gap_type_aware():
    sg = _irregular_signed()
    g_signed = sg.get_gap(laplacian_type=SG_LAPL_SIGNED)
    g_rw = sg.get_gap(laplacian_type=SG_LAPL_RW)
    g_sym = sg.get_gap(laplacian_type=SG_LAPL_SYM)
    assert g_rw == pytest.approx(g_sym, abs=1e-9)          # isospectral
    assert abs(g_signed - g_sym) > 1e-6                    # genuinely distinct
    # back to signed reproduces the original gap
    assert sg.get_gap(laplacian_type=SG_LAPL_SIGNED) == pytest.approx(g_signed, abs=1e-12)


# --------------------------------------------------------------------------- #
# 12. lrg free helper get_graph_lspectrum dispatches on laplacian_type
# --------------------------------------------------------------------------- #
def test_get_graph_lspectrum_dispatch():
    g = _weighted(nx.path_graph(7))
    g[0][1]["weight"] = -1.0
    _, w_sym = get_graph_lspectrum(g, laplacian_type=SG_LAPL_SYM)
    _, w_rw = get_graph_lspectrum(g, laplacian_type=SG_LAPL_RW)
    assert not np.iscomplexobj(w_rw)
    np.testing.assert_allclose(np.sort(w_sym), np.sort(w_rw), atol=1e-9)
    # back-compat alias returns the symmetric spectrum
    _, w_alias = get_graph_lspectrum_rw(g)
    np.testing.assert_allclose(np.sort(np.real(w_alias)), np.sort(w_sym), atol=1e-9)


# --------------------------------------------------------------------------- #
# 13. L2D pipeline kernel honours laplacian_type
# --------------------------------------------------------------------------- #
def test_l2d_kernel_laplacian_type():
    import sys
    from pathlib import Path

    src = Path(__file__).resolve().parents[2] / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    from kernels.L2D import eigv_for_lattice2D

    ev_signed = np.sort(eigv_for_lattice2D(
        side=6, pflip=0.1, geo="squared", laplacian_type=SG_LAPL_SIGNED))
    ev_rw = np.sort(eigv_for_lattice2D(
        side=6, pflip=0.1, geo="squared", laplacian_type=SG_LAPL_RW))
    # rw spectrum bounded in [0, 2]; signed is not (combinatorial)
    assert ev_rw.min() >= -1e-8 and ev_rw.max() <= 2.0 + 1e-8
    assert ev_signed.max() > 2.0  # combinatorial signed Laplacian of a degree-4 lattice


# --------------------------------------------------------------------------- #
# 14. graph-tool engine parity (skipped if graph-tool is unavailable)
# --------------------------------------------------------------------------- #
def test_gt_parity():
    pytest.importorskip("graph_tool")
    from lrgsglib.graphs.gt.SignedGraphGT.SignedGraphGT import SignedGraphGT

    g = _weighted(nx.path_graph(7))
    g[0][1]["weight"] = -1.0

    gt = SignedGraphGT.from_networkx(g)
    Lrw = gt.get_signed_rw_laplacian(sym=False)
    Lsym = gt.get_signed_rw_laplacian(sym=True)
    assert not np.allclose(Lrw, Lrw.T)
    assert np.allclose(Lsym, Lsym.T)

    gt.compute_laplacian_spectrum(laplacian_type=SG_LAPL_RW)
    ev_gt = np.sort(gt.eigv.copy())

    nxg = _signed_graph(g)
    nxg.compute_laplacian_spectrum(laplacian_type=SG_LAPL_RW)
    ev_nx = np.sort(nxg.eigv.copy())
    np.testing.assert_allclose(ev_gt, ev_nx, atol=1e-8)
