"""Signed heat kernel, synchronization partitions and planted-label helpers."""

import numpy as np
import pytest

from lrgsglib.graphs import StochasticBlockModel

SIZES = [30, 30, 30]
P_MATRIX = [[0.4, 0.1, 0.1], [0.1, 0.4, 0.1], [0.1, 0.1, 0.4]]
TAU = 5.0


@pytest.fixture(scope="module")
def sbm():
    g = StochasticBlockModel(sizes=SIZES, p_matrix=P_MATRIX, seed=7)
    g.flip_inter_block_edges()
    g.compute_laplacian_spectrum_weigV()
    return g


def test_planted_labels_and_block_signs(sbm):
    lab = sbm.get_planted_labels()
    assert lab.shape == (sbm.N,) and set(np.bincount(lab)) == {30}
    nodes = list(sbm.gr["G"].nodes())
    pos = {u: i for i, u in enumerate(nodes)}
    for u, v, w in sbm.gr["G"].edges(data="weight"):
        assert (w < 0) == (lab[pos[u]] != lab[pos[v]])
    assert len(sbm.fleset["G"]) == sbm.Ne_n


def test_kernel_identities(sbm):
    K = sbm.get_heat_kernel(TAU)
    Y = sbm.get_kernel_embedding(TAU)
    rho = sbm.get_kernel_correlation(TAU)
    assert np.allclose(K, Y @ Y.T)
    assert np.allclose(K, K.T)
    assert np.allclose(np.diag(rho), 1.0)
    assert rho.min() >= -1 - 1e-12 and rho.max() <= 1 + 1e-12
    # shift changes nothing but an overall factor
    K0 = sbm.get_heat_kernel(TAU, shift=False)
    assert np.allclose(K0 * np.exp(TAU * sbm.eigv[0]), K)


@pytest.mark.parametrize("kw", [dict(method="kmeans"), dict(method="linkage"),
                                dict(method="kmeans", mod_sign=True)])
def test_kernel_partition_recovers_blocks(sbm, kw):
    labels = sbm.get_kernel_partition(TAU, n_clusters=3, seed=0, **kw)
    assert labels.shape == (sbm.N,)
    assert sbm.score_partition(labels)["nmi"] > 0.9


def test_spectral_partition_and_coarse_graph(sbm):
    labels = sbm.get_spectral_partition(3, seed=0)
    assert sbm.score_partition(labels)["ari"] > 0.9
    # cached full spectrum untouched by the local bottom-k solve
    assert np.shape(sbm.eigV) == (sbm.N, sbm.N)
    Wp, Wm = sbm.get_coarse_signed_adjacency(sbm.get_planted_labels())
    assert Wp.shape == (3, 3) and np.all(np.diag(Wm) == 0)
    assert np.all(Wp[~np.eye(3, dtype=bool)] == 0)
    assert Wp.sum() / 2 + np.diag(Wp).sum() / 2 + Wm.sum() / 2 == sbm.Ne


def test_tau_peaks(sbm):
    peaks = sbm.get_tau_peaks()
    assert peaks.ndim == 1
    assert np.all(peaks > 0)
