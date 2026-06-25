"""Canonical voter / coarsening observables (Phase 1).

Unit tests for the engine-agnostic free functions in ``lrgsglib.utils.statsys``
(consensus/exit-time, survival, interface-density ensemble, susceptibility,
torus-wrap spanning) plus integration tests for the thin ``VoterModel`` ensemble
methods that drive them over real runs.
"""
import numpy as np
import pytest

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.VoterModel.VoterModel import VoterModel
from lrgsglib.utils.statsys import (
    cluster_wraps,
    consensus_time_stats,
    giant_cluster_spans,
    interface_density_ensemble,
    order_parameter_susceptibility,
    signed_neighbor_arrays,
    survival_curve,
)


# --------------------------------------------------------------------------- #
# Pure free functions (deterministic)
# --------------------------------------------------------------------------- #
@pytest.mark.code
def test_consensus_time_stats_basic():
    stats = consensus_time_stats([10, 20, None, 30], n_steps=100)
    assert stats["n"] == 4
    assert stats["censored"] == 1
    assert stats["p_consensus"] == pytest.approx(0.75)
    assert stats["mean"] == pytest.approx(20.0)
    assert stats["median"] == pytest.approx(20.0)
    assert stats["n_steps"] == 100


@pytest.mark.code
def test_consensus_time_stats_all_censored_is_nan():
    stats = consensus_time_stats([None, None, None], n_steps=50)
    assert stats["p_consensus"] == 0.0
    assert stats["censored"] == 3
    assert np.isnan(stats["mean"]) and np.isnan(stats["median"])


@pytest.mark.code
def test_survival_curve_monotone_and_S0():
    sc = survival_curve([2, 4, None], n_steps=6)
    assert sc[0] == pytest.approx(1.0)               # nobody absorbs before t=0
    assert np.all(np.diff(sc) <= 1e-12)              # monotone non-increasing
    # one run never absorbs -> tail floor = 1/3
    assert sc[-1] == pytest.approx(1.0 / 3.0)


@pytest.mark.code
def test_survival_curve_complements_consensus_fraction():
    times = [1, 3, 5, None]
    n_steps = 8
    sc = survival_curve(times, n_steps)
    stats = consensus_time_stats(times, n_steps)
    # by the horizon, the survivors are exactly the censored runs
    assert sc[-1] == pytest.approx(1.0 - stats["p_consensus"])


@pytest.mark.code
def test_interface_density_ensemble_forward_fills():
    mean, std = interface_density_ensemble([[0.5, 0.2, 0.0], [0.4, 0.1]])
    # second run is forward-filled at 0.1 -> mean = [0.45, 0.15, 0.05]
    assert mean == pytest.approx([0.45, 0.15, 0.05])
    assert std.shape == (3,)


@pytest.mark.code
def test_interface_density_ensemble_single_run_is_identity():
    series = [0.6, 0.4, 0.3, 0.1]
    mean, std = interface_density_ensemble([series])
    assert mean == pytest.approx(series)
    assert np.allclose(std, 0.0)


@pytest.mark.code
def test_order_parameter_susceptibility_matches_N_var():
    samples = [0.9, 0.1, 0.5, 0.5]
    mean, chi = order_parameter_susceptibility(samples, n_sites=100)
    assert mean == pytest.approx(0.5)
    assert chi == pytest.approx(100 * np.var(samples))   # population variance


@pytest.mark.code
def test_cluster_wraps_row_spans_blob_does_not():
    L = 6
    row = np.zeros((L, L), bool)
    row[0, :] = True                                  # a full row wraps x
    blob = np.zeros((L, L), bool)
    blob[1:3, 1:3] = True                             # a 2x2 blob does not
    assert cluster_wraps(row) is True
    assert cluster_wraps(blob) is False
    assert cluster_wraps(np.zeros((L, L), bool)) is False


@pytest.mark.code
def test_cluster_wraps_3d():
    L = 4
    col = np.zeros((L, L, L), bool)
    col[:, 0, 0] = True                               # a full line wraps axis 0
    assert cluster_wraps(col) is True


# --------------------------------------------------------------------------- #
# giant_cluster_spans on a real lattice substrate
# --------------------------------------------------------------------------- #
@pytest.mark.code
def test_giant_cluster_spans_consensus_vs_checkerboard():
    lat = Lattice2D(8, engine="nx")
    idx, b = signed_neighbor_arrays(lat, "rawspin")
    consensus = np.ones(lat.N, dtype=np.int8)
    assert giant_cluster_spans(consensus, idx, b, lat.syshape) is True
    # checkerboard: no two equal-spin neighbours -> every node its own cluster
    rows, cols = lat.syshape
    chk = np.fromfunction(lambda i, j: ((i + j) % 2) * 2 - 1, (rows, cols))
    chk = chk.astype(np.int8).ravel()
    assert giant_cluster_spans(chk, idx, b, lat.syshape) is False


# --------------------------------------------------------------------------- #
# VoterModel ensemble methods over real runs
# --------------------------------------------------------------------------- #
def _run_ensemble(lat_factory, n_runs, **vm_kw):
    runs = []
    for seed in range(n_runs):
        lat = lat_factory(seed)
        vm = VoterModel(sg=lat, seed=seed, savedisk=False, runlang="py", **vm_kw)
        vm.run()
        runs.append(vm)
    return runs


@pytest.mark.physical
def test_voter_consensus_time_unsigned_absorbs():
    # Unsigned connected lattice: an absorbing (consensus) state exists, so with a
    # generous horizon every run freezes -> p_consensus == 1, finite mean time.
    runs = _run_ensemble(
        lambda s: Lattice2D(6, engine="nx"), n_runs=5,
        steps=20000, absorbing_check=True, ic="random",
    )
    stats = VoterModel.consensus_time_stats(runs)
    assert stats["n"] == 5
    assert stats["p_consensus"] == pytest.approx(1.0)
    assert stats["mean"] > 0.0
    # survival curve is a proper monotone S(t) ending at 0 (all absorbed)
    sc = VoterModel.survival_curve(runs)
    assert sc[0] == pytest.approx(1.0)
    assert np.all(np.diff(sc) <= 1e-12)


@pytest.mark.physical
def test_voter_interface_density_ensemble_needs_trajectory():
    runs = _run_ensemble(
        lambda s: Lattice2D(6, engine="nx"), n_runs=3,
        steps=40, savedyn=True, ic="random",
    )
    mean, std = VoterModel.interface_density_ensemble(runs)
    assert mean.ndim == 1 and mean.size > 0
    assert np.all(mean >= -1e-9) and np.all(mean <= 1.0 + 1e-9)


@pytest.mark.physical
def test_voter_order_parameter_susceptibility_wires_up():
    runs = _run_ensemble(
        lambda s: Lattice2D(6, engine="nx"), n_runs=4,
        steps=200, ic="random",
    )
    mean, chi = VoterModel.order_parameter_susceptibility(runs)
    assert 0.0 <= mean <= 1.0
    assert chi >= 0.0


@pytest.mark.code
def test_voter_giant_cluster_spans_method():
    lat = Lattice2D(8, engine="nx")
    vm = VoterModel(sg=lat, steps=10, savedisk=False, runlang="py")
    assert vm.giant_cluster_spans(s=np.ones(lat.N, dtype=np.int8)) is True
