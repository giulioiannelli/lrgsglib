"""Cluster-size-distribution observable for the VoterModel.

A cluster is a connected component of the subgraph of *active* edges:
``b_ij * s_i * s_j > 0`` with ``b_ij = sign(w_ij)`` ('satisfied' mode) or
``b_ij = +1`` ('rawspin' mode). The incremental ``ClusterTracker`` (union-find
merges + bounded local rescan on splits) is validated against the brute-force
``cluster_size_distribution`` oracle, both standalone (random flip sequences) and
wired into the Python single-spin samplers (asynchronous / gillespie).

Signed-substrate / absorbing physics references: ``iannelli2025topological.pdf``
(signed Laplacian, frustration); cluster semantics in ``VoterModel/defaults.py``.
"""

import numpy as np
import pytest

from lrgsglib.statsys.VoterModel._cluster_tracker import (
    ClusterTracker,
    cluster_size_distribution,
    edge_sign_arrays,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _random_signed_graph(n, p, neg_frac, seed):
    """Ragged (idx, signs) for a random signed graph (signs = sign(w_ij))."""
    rng = np.random.default_rng(seed)
    idx = [[] for _ in range(n)]
    signs = [[] for _ in range(n)]
    for u in range(n):
        for v in range(u + 1, n):
            if rng.random() < p:
                w = -1 if rng.random() < neg_frac else 1
                idx[u].append(v); signs[u].append(w)
                idx[v].append(u); signs[v].append(w)
    idx = [np.asarray(a, dtype=np.int64) for a in idx]
    signs = [np.asarray(a, dtype=np.int8) for a in signs]
    return idx, signs


def _dist(s, idx, b):
    return {int(k): int(v) for k, v in cluster_size_distribution(s, idx, b).items()}


def _lat(tmp_path, side=12, pflip=0.0, seed=1):
    from lrgsglib.graphs.nx import Lattice2DNX
    return Lattice2DNX(
        side1=side, geo="sqr", pbc=True, pflip=pflip, seed=seed,
        path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
    )


def _model_oracle(vm):
    idx, signs = vm._gillespie_neighbors()
    b = edge_sign_arrays(signs, vm.cluster_mode)
    return _dist(vm.s, idx, b)


# ===================================================================
# Incremental tracker vs brute-force oracle (random flip sequences)
# ===================================================================
@pytest.mark.parametrize("mode", ["satisfied", "rawspin"])
@pytest.mark.parametrize("neg_frac", [0.0, 0.4, 0.7])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_tracker_matches_oracle_random(mode, neg_frac, seed):
    n, p, n_flips = 40, 0.18, 120
    idx, signs = _random_signed_graph(n, p, neg_frac, seed)
    b = edge_sign_arrays(signs, mode)
    rng = np.random.default_rng(seed + 7)
    s = rng.choice(np.array([-1, 1], dtype=np.int8), size=n).astype(np.int8)
    tr = ClusterTracker(s, idx, b)
    assert tr.size_distribution() == _dist(s, idx, b)  # initial partition
    for _ in range(n_flips):
        i = int(rng.integers(n))
        tr.flip(i)                          # pre-flip bookkeeping
        s[i] = np.int8(-int(s[i]))          # the model's flip
        assert tr.size_distribution() == _dist(s, idx, b)
        assert sum(sz * c for sz, c in tr.size_distribution().items()) == n


def test_tracker_degenerate_cases():
    """Singleton removal, leaf removal (no split), and a real split."""
    # Path 0-1-2-3-4 (all +1, all positive edges) is one cluster of size 5.
    n = 5
    idx = [np.array(a, dtype=np.int64) for a in
           ([1], [0, 2], [1, 3], [2, 4], [3])]
    signs = [np.ones(len(a), dtype=np.int8) for a in idx]
    b = edge_sign_arrays(signs, "satisfied")
    s = np.ones(n, dtype=np.int8)
    tr = ClusterTracker(s, idx, b)
    assert tr.size_distribution() == {5: 1}
    # Flip the middle node 2 (it has 2 agreeing neighbours -> path splits 0-1 | 3-4
    # and node 2 becomes its own singleton).
    tr.flip(2); s[2] = np.int8(-1)
    assert tr.size_distribution() == _dist(s, idx, b)
    assert tr.size_distribution() == {2: 2, 1: 1}
    # Flip a leaf (node 0): no split, it just leaves its 2-cluster.
    tr.flip(0); s[0] = np.int8(-1)
    assert tr.size_distribution() == _dist(s, idx, b)


# ===================================================================
# Wired into the Python single-spin samplers
# ===================================================================
@pytest.mark.physical
@pytest.mark.parametrize("upd_mode", ["asynchronous", "gillespie"])
@pytest.mark.parametrize("cluster_mode", ["satisfied", "rawspin"])
@pytest.mark.parametrize("pflip", [0.0, 0.2])
def test_model_cluster_dist_tracks(tmp_path, upd_mode, cluster_mode, pflip):
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=14, pflip=pflip, seed=5)
    lat.flip_random_fract_edges()
    vm = VoterModel(
        sg=lat, steps=30, runlang="py", rule="linear", upd_mode=upd_mode,
        track_clusters=True, cluster_mode=cluster_mode,
        save_magnetization=True, absorbing_check=(pflip == 0.0), seed=3,
    )
    vm.run(tqdm_on=False)
    # one distribution recorded per recorded magnetization sample
    assert len(vm.cluster_dist) == len(vm.magn)
    # every recorded distribution conserves the node count
    for dist in vm.cluster_dist:
        assert sum(sz * c for sz, c in dist.items()) == vm.N
    # the live tracker is exactly in sync with an independent oracle
    assert vm._tracker.size_distribution() == _model_oracle(vm)


@pytest.mark.physical
def test_model_cluster_dist_nonlinear_async(tmp_path):
    """Cluster tracking works for any rule under the asynchronous schedule."""
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=12, pflip=0.15, seed=9)
    lat.flip_random_fract_edges()
    vm = VoterModel(
        sg=lat, steps=25, runlang="py", rule="majority",
        upd_mode="asynchronous", track_clusters=True, seed=1,
    )
    vm.run(tqdm_on=False)
    assert vm._tracker.size_distribution() == _model_oracle(vm)


# ===================================================================
# Guards: backends / schedules that do not (yet) emit clusters
# ===================================================================
def test_track_clusters_rejects_sync_and_link(tmp_path):
    from lrgsglib.statsys import VoterModel

    for upd in ("synchronous", "link"):
        vm = VoterModel(sg=_lat(tmp_path), steps=5, runlang="py",
                        upd_mode=upd, track_clusters=True, seed=1)
        with pytest.raises(NotImplementedError):
            vm.run(tqdm_on=False)


def test_track_clusters_rejects_unsupported_native(tmp_path):
    """Native cluster tracking is pybind+gillespie only; other native combos raise."""
    from lrgsglib.statsys import VoterModel

    # pybind with a non-gillespie schedule
    vm = VoterModel(sg=_lat(tmp_path), steps=5, runlang="pb_voter",
                    upd_mode="asynchronous", track_clusters=True, seed=1)
    with pytest.raises(NotImplementedError):
        vm.run(tqdm_on=False)
    # C subprocess does not emit clusters yet
    vm = VoterModel(sg=_lat(tmp_path), steps=5, runlang="C0",
                    upd_mode="gillespie", track_clusters=True, seed=1)
    with pytest.raises(NotImplementedError):
        vm.run(tqdm_on=False)


@pytest.mark.physical
@pytest.mark.parametrize("cluster_mode", ["satisfied", "rawspin"])
@pytest.mark.parametrize("pflip", [0.0, 0.2])
def test_native_pb_gillespie_clusters(tmp_path, cluster_mode, pflip):
    """pybind gillespie emits a cluster distribution matching the oracle."""
    pytest.importorskip("lrgsglib.statsys.VoterModel.ccore", reason="pb build")
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=16, pflip=pflip, seed=11)
    lat.flip_random_fract_edges()
    vm = VoterModel(
        sg=lat, steps=40, runlang="pb_voter", rule="linear", upd_mode="gillespie",
        track_clusters=True, cluster_mode=cluster_mode,
        save_magnetization=True, absorbing_check=(pflip == 0.0), seed=4,
    )
    try:
        vm.run(tqdm_on=False)
    except (ImportError, RuntimeError, OSError):
        pytest.skip("native voter module unavailable")
    assert len(vm.cluster_dist) == len(vm.magn)
    for dist in vm.cluster_dist:
        assert sum(sz * c for sz, c in dist.items()) == vm.N
    # final recorded distribution equals an independent oracle on the final state
    assert vm.cluster_dist[-1] == _model_oracle(vm)


def test_invalid_cluster_mode_raises(tmp_path):
    from lrgsglib.statsys import VoterModel

    with pytest.raises(ValueError):
        VoterModel(sg=_lat(tmp_path), steps=5, runlang="py",
                   track_clusters=True, cluster_mode="bogus", seed=1)
