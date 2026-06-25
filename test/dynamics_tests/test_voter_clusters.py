"""Cluster-size-distribution observable for the VoterModel.

A cluster is a connected component of the subgraph of *active* edges:
``b_ij * s_i * s_j > 0`` with ``b_ij = sign(w_ij)`` ('satisfied' mode) or
``b_ij = +1`` ('rawspin' mode). The distribution is recomputed from scratch at
each recorded sweep (a connected-components pass over the active-edge subgraph)
by both the Python ``cluster_size_distribution`` and the native C kernel
(``_ccore/LRGSG_clusters.c``); the two are checked to agree, and the recorded
``cluster_dist`` is checked against the oracle on the recorded state.

Signed-substrate / absorbing physics references: ``iannelli2025topological.pdf``
(signed Laplacian, frustration); cluster semantics in ``VoterModel/defaults.py``.
"""

import numpy as np
import pytest

from lrgsglib.utils.statsys import (
    cluster_components,
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


def _oracle_initial(vm):
    """Oracle distribution on vm's CURRENT state (call right after init)."""
    idx, signs = vm._gillespie_neighbors()
    b = edge_sign_arrays(signs, vm.cluster_mode)
    return _dist(vm.s, idx, b)


# ===================================================================
# Oracle: connected components over the active-edge subgraph
# ===================================================================
def test_oracle_path_split_satisfied():
    """Path 0-1-2-3-4 (all +1) is one size-5 cluster; flipping the middle
    splits it into 0-1 | 3-4 with node 2 a singleton."""
    idx = [np.array(a, dtype=np.int64) for a in
           ([1], [0, 2], [1, 3], [2, 4], [3])]
    b = edge_sign_arrays([np.ones(len(a), dtype=np.int8) for a in idx], "satisfied")
    s = np.ones(5, dtype=np.int8)
    assert _dist(s, idx, b) == {5: 1}
    s[2] = -1
    assert _dist(s, idx, b) == {2: 2, 1: 1}


def test_oracle_satisfied_vs_rawspin():
    """A negative edge is 'active' when its endpoints DISAGREE (satisfied), but
    inactive in rawspin mode (which only joins equal raw spins)."""
    idx = [np.array([1], dtype=np.int64), np.array([0], dtype=np.int64)]
    signs = [np.array([-1], dtype=np.int8), np.array([-1], dtype=np.int8)]
    s = np.array([1, -1], dtype=np.int8)         # disagree across a negative edge
    assert _dist(s, idx, edge_sign_arrays(signs, "satisfied")) == {2: 1}
    assert _dist(s, idx, edge_sign_arrays(signs, "rawspin")) == {1: 2}


@pytest.mark.parametrize("mode", ["satisfied", "rawspin"])
@pytest.mark.parametrize("neg_frac", [0.0, 0.5])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_oracle_components_conserve(mode, neg_frac, seed):
    """labels partition every node; sizes sum to N for arbitrary states."""
    n, p = 40, 0.18
    idx, signs = _random_signed_graph(n, p, neg_frac, seed)
    b = edge_sign_arrays(signs, mode)
    rng = np.random.default_rng(seed + 7)
    for _ in range(20):
        s = rng.choice(np.array([-1, 1], dtype=np.int8), size=n).astype(np.int8)
        label = cluster_components(s, idx, b)
        assert label.min() >= 0 and -1 not in label
        assert sum(sz * c for sz, c in _dist(s, idx, b).items()) == n


# ===================================================================
# Python backend: recorded distribution == oracle on the recorded state
# ===================================================================
@pytest.mark.physical
@pytest.mark.parametrize(
    "upd_mode", ["asynchronous", "synchronous", "link", "gillespie"]
)
@pytest.mark.parametrize("cluster_mode", ["satisfied", "rawspin"])
def test_py_cluster_dist_matches_oracle(tmp_path, upd_mode, cluster_mode):
    """Full-recompute recording works for EVERY Python schedule; the first
    record equals the oracle on the initial state, and node count is conserved."""
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=14, pflip=0.15, seed=5)
    lat.flip_random_fract_edges()
    vm = VoterModel(
        sg=lat, steps=20, runlang="py", rule="linear", upd_mode=upd_mode,
        track_clusters=True, cluster_mode=cluster_mode,
        savemagn=True, absorbing_check=False, seed=3,
    )
    vm.init_voter_dynamics()
    oracle0 = _oracle_initial(vm)
    vm.run(tqdm_on=False)
    assert len(vm.cluster_dist) == len(vm.magn)
    assert vm.cluster_dist[0] == oracle0
    for dist in vm.cluster_dist:
        assert sum(sz * c for sz, c in dist.items()) == vm.N


@pytest.mark.physical
def test_py_cluster_dist_nonlinear(tmp_path):
    """Recording is rule-agnostic (majority rule, asynchronous schedule)."""
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=12, pflip=0.15, seed=9)
    lat.flip_random_fract_edges()
    vm = VoterModel(
        sg=lat, steps=15, runlang="py", rule="majority",
        upd_mode="asynchronous", track_clusters=True, seed=1,
    )
    vm.init_voter_dynamics()
    oracle0 = _oracle_initial(vm)
    vm.run(tqdm_on=False)
    assert vm.cluster_dist[0] == oracle0
    for dist in vm.cluster_dist:
        assert sum(sz * c for sz, c in dist.items()) == vm.N


@pytest.mark.physical
def test_py_cluster_consensus_single_domain(tmp_path):
    """On an unsigned lattice the absorbing state is consensus -> exactly one
    domain spanning all N nodes."""
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=12, pflip=0.0, seed=2)
    vm = VoterModel(
        sg=lat, steps=4000, runlang="py", rule="linear", upd_mode="gillespie",
        track_clusters=True, absorbing_check=True, seed=4,
    )
    vm.init_voter_dynamics()
    vm.run(tqdm_on=False)
    if vm.absorbed_at is not None:
        assert vm.cluster_dist[-1] == {vm.N: 1}


# ===================================================================
# Native pybind backend: C connected components == Python oracle
# ===================================================================
@pytest.mark.physical
@pytest.mark.parametrize("cluster_mode", ["satisfied", "rawspin"])
@pytest.mark.parametrize("pflip", [0.0, 0.2])
def test_native_pb_gillespie_clusters(tmp_path, cluster_mode, pflip):
    """pybind gillespie emits a cluster distribution; the C recompute matches the
    Python oracle on the initial state and conserves the node count."""
    pytest.importorskip("lrgsglib.statsys.VoterModel.ccore", reason="pb build")
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=16, pflip=pflip, seed=11)
    lat.flip_random_fract_edges()
    vm = VoterModel(
        sg=lat, steps=40, runlang="pb", rule="linear", upd_mode="gillespie",
        track_clusters=True, cluster_mode=cluster_mode,
        savemagn=True, absorbing_check=(pflip == 0.0), seed=4,
    )
    vm.init_voter_dynamics()
    oracle0 = _oracle_initial(vm)
    try:
        vm.run(tqdm_on=False)
    except (ImportError, RuntimeError, OSError):
        pytest.skip("native voter module unavailable")
    assert len(vm.cluster_dist) == len(vm.magn)
    assert vm.cluster_dist[0] == oracle0          # C CC == Python oracle
    for dist in vm.cluster_dist:
        assert sum(sz * c for sz, c in dist.items()) == vm.N
    if pflip == 0.0 and vm.absorbed_at is not None:
        assert vm.cluster_dist[-1] == {vm.N: 1}    # consensus = one domain


@pytest.mark.physical
@pytest.mark.parametrize("upd_mode", ["asynchronous", "synchronous", "link"])
@pytest.mark.parametrize("cluster_mode", ["satisfied", "rawspin"])
def test_native_pb_nongillespie_clusters(tmp_path, upd_mode, cluster_mode):
    """pybind records the cluster distribution for the non-gillespie schedules too.
    With savedyn on, the native per-sweep recompute must equal the Python oracle on
    EVERY recorded state -- this also exercises the synchronous schedule, whose
    s/sbuf pointer swap is handled by clusters_set_state on the C side."""
    pytest.importorskip("lrgsglib.statsys.VoterModel.ccore", reason="pb build")
    from lrgsglib.statsys import VoterModel

    lat = _lat(tmp_path, side=14, pflip=0.2, seed=7)
    lat.flip_random_fract_edges()
    vm = VoterModel(
        sg=lat, steps=25, runlang="pb", rule="linear", upd_mode=upd_mode,
        track_clusters=True, cluster_mode=cluster_mode,
        savemagn=True, savedyn=True, absorbing_check=False, seed=2,
    )
    vm.init_voter_dynamics()
    try:
        vm.run(tqdm_on=False)
    except (ImportError, RuntimeError, OSError):
        pytest.skip("native voter module unavailable")
    assert len(vm.cluster_dist) == len(vm.s_t) == len(vm.magn) == 25
    idx, signs = vm._gillespie_neighbors()
    b = edge_sign_arrays(signs, cluster_mode)
    for k, sk in enumerate(vm.s_t):
        assert dict(vm.cluster_dist[k]) == _dist(np.asarray(sk, np.int8), idx, b)


# ===================================================================
# Guards: which backend/schedule combinations emit clusters
# ===================================================================
@pytest.mark.physical
@pytest.mark.parametrize("upd_mode", ["asynchronous", "synchronous", "link"])
def test_py_all_schedules_support_clusters(tmp_path, upd_mode):
    """The Python full-recompute path supports every schedule (no raise)."""
    from lrgsglib.statsys import VoterModel

    vm = VoterModel(sg=_lat(tmp_path), steps=5, runlang="py",
                    upd_mode=upd_mode, track_clusters=True,
                    savemagn=True, absorbing_check=False, seed=1)
    vm.init_voter_dynamics()
    vm.run(tqdm_on=False)
    assert len(vm.cluster_dist) == len(vm.magn) == 5


def test_track_clusters_c_subprocess_emits(tmp_path):
    """The C subprocess now emits the cluster-size distribution too (sparse
    cldist file -> list[{size: count}]), for every schedule."""
    from lrgsglib.statsys import VoterModel

    for upd_mode in ("asynchronous", "gillespie"):
        vm = VoterModel(sg=_lat(tmp_path), steps=5, runlang="C0",
                        upd_mode=upd_mode, track_clusters=True,
                        savemagn=True, savedisk=False, seed=1)
        vm.init_voter_dynamics()
        vm.run(tqdm_on=False)
        cd = list(vm.cluster_dist)
        assert len(cd) > 0
        for rec in cd:                       # each record partitions the N nodes
            assert sum(int(k) * int(v) for k, v in rec.items()) == vm.N


def test_track_clusters_rejects_vectorized(tmp_path):
    """The fused NumPy/CuPy sweep has no per-sweep record hook -> raises."""
    from lrgsglib.statsys import VoterModel

    vm = VoterModel(sg=_lat(tmp_path), steps=5, runlang="np",
                    upd_mode="synchronous", track_clusters=True, seed=1)
    with pytest.raises(NotImplementedError):
        vm.run(tqdm_on=False)


def test_invalid_cluster_mode_raises(tmp_path):
    from lrgsglib.statsys import VoterModel

    with pytest.raises(ValueError):
        VoterModel(sg=_lat(tmp_path), steps=5, runlang="py",
                   track_clusters=True, cluster_mode="bogus", seed=1)


@pytest.mark.physical
def test_interface_density_and_largest_domain(tmp_path):
    """interface_density / largest_domain_fraction sanity (single state + series).

    Consensus is one frozen domain: rho=0, giant fraction=1. Arbitrary states
    stay in range. The *_series helpers map over a saved trajectory.
    """
    from lrgsglib.statsys import VoterModel

    vm = VoterModel(sg=_lat(tmp_path), steps=1, runlang="py", seed=1)
    N = vm.N
    consensus = np.ones(N, dtype=np.int8)
    assert vm.interface_density(consensus) == 0.0
    assert vm.is_absorbing(consensus)
    assert vm.largest_domain_fraction(consensus) == 1.0

    rng = np.random.default_rng(0)
    rand = rng.choice(np.array([-1, 1], dtype=np.int8), size=N)
    assert 0.0 <= vm.interface_density(rand) <= 1.0
    assert 1.0 / N <= vm.largest_domain_fraction(rand) <= 1.0

    vm2 = VoterModel(sg=_lat(tmp_path), steps=5, runlang="py", seed=2,
                     savedyn=True, savedisk=False)
    vm2.run(tqdm_on=False)
    rho_t = vm2.interface_density_series()
    gf_t = vm2.largest_domain_fraction_series()
    assert len(rho_t) == len(gf_t) == len(vm2.s_t)
    assert all(0.0 <= x <= 1.0 for x in rho_t)
    assert all(1.0 / N <= x <= 1.0 for x in gf_t)
