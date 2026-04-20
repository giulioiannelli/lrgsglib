"""Statistical parity: Python ``run_walker`` vs pybind11 ``pb_absorb``.

Both backends use their own RNG streams (numpy PCG64 for Python, SFMT
for C) so per-walker outputs differ; the ensemble-averaged observables
must agree within a few SEM.  A visits-count cosine similarity > 0.9
confirms the aggregated visit field is in the same direction.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.SignedRW import AbsorbingWalker
from lrgsglib.statsys.SignedRW.SignedWalker import _load_srw_native

pytestmark = pytest.mark.skipif(
    _load_srw_native() is None,
    reason="_srw_native extension not built; run make under ccore/native/",
)


@pytest.fixture(scope="module")
def lat_dis():
    tmp = tempfile.TemporaryDirectory(prefix="srw_parity_")
    lat = Lattice2D(32, geo='sqr', pbc=True, pflip=0.10, seed=1,
                    path_data=Path(tmp.name))
    lat.flip_random_fract_edges()
    yield lat
    tmp.cleanup()


def _run(lat, runlang, seed, n_walkers=5000, coverage_stop=0.05):
    w = AbsorbingWalker(lat, n_walkers=n_walkers, seed=seed,
                        coverage_stop=coverage_stop, runlang=runlang)
    w.run()
    return w


@pytest.mark.parametrize("seed", [1, 7, 42])
def test_mean_stop_step_parity(lat_dis, seed):
    w_py = _run(lat_dis, 'py', seed)
    w_pb = _run(lat_dis, 'pb_absorb', seed)
    # Compare means: |Δmean| must be < 4 SEM summed in quadrature.
    mean_py = w_py.stop_step.mean()
    mean_pb = w_pb.stop_step.mean()
    sem = np.hypot(w_py.stop_step.std(ddof=1) / np.sqrt(len(w_py.stop_step)),
                   w_pb.stop_step.std(ddof=1) / np.sqrt(len(w_pb.stop_step)))
    assert abs(mean_py - mean_pb) < 4 * sem, (
        f"mean_step py={mean_py:.3f} pb={mean_pb:.3f} diff/SEM="
        f"{abs(mean_py - mean_pb) / sem:.2f}"
    )


@pytest.mark.parametrize("seed", [1, 7, 42])
def test_mean_unique_frac_parity(lat_dis, seed):
    w_py = _run(lat_dis, 'py', seed)
    w_pb = _run(lat_dis, 'pb_absorb', seed)
    mean_py = w_py.unique_frac.mean()
    mean_pb = w_pb.unique_frac.mean()
    # Absolute tolerance 5% of the py value (generous — 5000 walkers have
    # low SEM but we allow some systematic drift from different RNGs).
    assert abs(mean_py - mean_pb) < 0.05 * mean_py + 1e-4


@pytest.mark.parametrize("seed", [1, 7, 42])
def test_visits_agg_cosine_parity(lat_dis, seed):
    w_py = _run(lat_dis, 'py', seed)
    w_pb = _run(lat_dis, 'pb_absorb', seed)
    c1 = w_py.visits_agg.astype(np.float64)
    c2 = w_pb.visits_agg.astype(np.float64)
    cos = c1.dot(c2) / (np.linalg.norm(c1) * np.linalg.norm(c2))
    assert cos > 0.9, f"visits_agg cosine similarity {cos:.4f} < 0.9"


@pytest.mark.parametrize("seed", [1, 7])
def test_death_count_parity(lat_dis, seed):
    """At low coverage on a disordered lattice, almost all walkers die."""
    w_py = _run(lat_dis, 'py', seed)
    w_pb = _run(lat_dis, 'pb_absorb', seed)
    died_py = int((w_py.stop_reason == 1).sum())
    died_pb = int((w_pb.stop_reason == 1).sum())
    # both should be > 99% of walkers
    n = len(w_py.stop_step)
    assert died_py / n > 0.99
    assert died_pb / n > 0.99


def test_runlang_rejects_non_absorb_rule(lat_dis):
    from lrgsglib.statsys.SignedRW import KillingWalker
    w = KillingWalker(lat_dis, n_walkers=10, seed=1,
                      coverage_stop=0.2, runlang='pb_absorb')
    with pytest.raises(RuntimeError, match="rule='absorb'"):
        w.run()


def test_runlang_rejects_unknown_pb(lat_dis):
    w = AbsorbingWalker(lat_dis, n_walkers=10, seed=1,
                        coverage_stop=0.2, runlang='pb_unknown')
    with pytest.raises(NotImplementedError, match="Phase 2"):
        w.run()


def test_seed_reproducibility_pb(lat_dis):
    """Same seed → identical output on the C backend."""
    w1 = _run(lat_dis, 'pb_absorb', seed=42)
    w2 = _run(lat_dis, 'pb_absorb', seed=42)
    np.testing.assert_array_equal(w1.stop_step, w2.stop_step)
    np.testing.assert_array_equal(w1.stop_reason, w2.stop_reason)
    np.testing.assert_array_equal(w1.visits_agg, w2.visits_agg)
    np.testing.assert_array_equal(w1.final_position, w2.final_position)
