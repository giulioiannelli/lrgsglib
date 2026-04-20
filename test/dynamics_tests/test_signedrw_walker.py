"""Unit tests for the :mod:`lrgsglib.statsys.SignedRW` walker family."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.SignedRW import (
    AbsorbingWalker,
    KillingWalker,
    SignedWalker,
    StickyWalker,
)
from lrgsglib.statsys.SignedRW._kernel import (
    run_walker,
    signed_lattice_tables,
)


@pytest.fixture(scope="module")
def tmp_data_dir():
    td = tempfile.TemporaryDirectory(prefix="srw_test_")
    yield Path(td.name)
    td.cleanup()


@pytest.fixture
def clean_lattice(tmp_data_dir):
    """Clean 16x16 sqr lattice with no negative edges (pflip=0)."""
    return Lattice2D(16, geo='sqr', pbc=True, pflip=0.0, seed=1,
                     path_data=tmp_data_dir)


@pytest.fixture
def disordered_lattice(tmp_data_dir):
    """16x16 sqr lattice with 10% negative edges."""
    lat = Lattice2D(16, geo='sqr', pbc=True, pflip=0.1, seed=1,
                    path_data=tmp_data_dir)
    lat.flip_random_fract_edges()
    return lat


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------

def test_signed_walker_is_abstract(clean_lattice):
    with pytest.raises(TypeError, match="abstract"):
        SignedWalker(clean_lattice)


# ---------------------------------------------------------------------------
# Smoke: each rule instantiates and runs end-to-end
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("Cls,rule", [
    (AbsorbingWalker, 'absorb'),
    (KillingWalker, 'kill'),
    (StickyWalker, 'sticky'),
])
def test_walker_runs(disordered_lattice, Cls, rule):
    w = Cls(disordered_lattice, n_walkers=30, seed=11,
            coverage_stop=0.2, x_node_behavior='reflect')
    w.run()
    assert w.rule == rule
    assert w.unique_visits.shape == (30,)
    assert w.stop_step.shape == (30,)
    assert w.stop_reason.shape == (30,)
    assert (w.unique_visits >= 1).all()
    assert (w.unique_visits <= disordered_lattice.N).all()
    assert set(np.unique(w.stop_reason)).issubset({0, 1, 2})
    assert w.visits_agg.sum() > 0
    # Every walker must have stopped.
    assert (w.stop_step >= 0).all()


# ---------------------------------------------------------------------------
# Seed reproducibility
# ---------------------------------------------------------------------------

def test_seed_reproducibility(disordered_lattice):
    w1 = AbsorbingWalker(disordered_lattice, n_walkers=50, seed=42,
                         coverage_stop=0.2)
    w1.run()
    w2 = AbsorbingWalker(disordered_lattice, n_walkers=50, seed=42,
                         coverage_stop=0.2)
    w2.run()
    np.testing.assert_array_equal(w1.stop_step, w2.stop_step)
    np.testing.assert_array_equal(w1.stop_reason, w2.stop_reason)
    np.testing.assert_array_equal(w1.final_position, w2.final_position)
    np.testing.assert_array_equal(w1.visits_agg, w2.visits_agg)


def test_different_seeds_differ(disordered_lattice):
    w1 = AbsorbingWalker(disordered_lattice, n_walkers=50, seed=1,
                         coverage_stop=0.2)
    w1.run()
    w2 = AbsorbingWalker(disordered_lattice, n_walkers=50, seed=2,
                         coverage_stop=0.2)
    w2.run()
    assert not np.array_equal(w1.stop_step, w2.stop_step)


# ---------------------------------------------------------------------------
# Start protocols
# ---------------------------------------------------------------------------

def test_start_fixed_single_node(disordered_lattice):
    w = AbsorbingWalker(disordered_lattice, n_walkers=20, seed=1,
                        start='fixed', start_kwargs={'node': 7},
                        coverage_stop=0.2)
    w.run()
    # every walker starts at node 7, so visits_agg[7] >= n_walkers
    assert w.visits_agg[7] >= w.n_walkers


def test_start_fixed_per_walker_array(disordered_lattice):
    nodes = np.arange(20, dtype=np.int64)
    w = AbsorbingWalker(disordered_lattice, n_walkers=20, seed=1,
                        start='fixed', start_kwargs={'nodes': nodes},
                        coverage_stop=0.2)
    w.run()
    assert (w.visits_agg[nodes] >= 1).all()


def test_start_center_2d(disordered_lattice):
    N = disordered_lattice.N
    side = int(round(N ** 0.5))
    expected = (side // 2) * (side + 1)
    w = AbsorbingWalker(disordered_lattice, n_walkers=10, seed=1,
                        start='center', coverage_stop=0.2)
    w.run()
    assert w.visits_agg[expected] >= w.n_walkers


def test_start_fixed_requires_node(disordered_lattice):
    with pytest.raises(ValueError, match="start_kwargs"):
        AbsorbingWalker(disordered_lattice, n_walkers=5, seed=1,
                        start='fixed', coverage_stop=0.2).run()


def test_bad_start_protocol(disordered_lattice):
    with pytest.raises(ValueError, match="start"):
        AbsorbingWalker(disordered_lattice, n_walkers=5, seed=1,
                        start='quantum_entangled', coverage_stop=0.2)


# ---------------------------------------------------------------------------
# X-node behaviour
# ---------------------------------------------------------------------------

def test_x_node_absorb_vs_reflect_differ(disordered_lattice):
    """'absorb' kills walkers faster than 'reflect' on a disordered lattice."""
    w_abs = AbsorbingWalker(disordered_lattice, n_walkers=300, seed=3,
                            coverage_stop=0.2, x_node_behavior='absorb')
    w_abs.run()
    w_ref = AbsorbingWalker(disordered_lattice, n_walkers=300, seed=3,
                            coverage_stop=0.2, x_node_behavior='reflect')
    w_ref.run()
    # On a non-trivial disorder, reflect-mode walkers live at least as long
    # as absorb-mode walkers (they can't die on unfrustrated negatives).
    assert w_ref.stop_step.mean() >= w_abs.stop_step.mean()


def test_clean_lattice_no_deaths(clean_lattice):
    """Zero disorder → no negative edges → no deaths; walkers end on coverage."""
    N = clean_lattice.N
    coverage = 0.2
    stop_thresh = max(1, int(coverage * N))
    w = AbsorbingWalker(clean_lattice, n_walkers=30, seed=1,
                        coverage_stop=coverage)
    w.run()
    assert (w.stop_reason == 0).all()
    # coverage threshold is floor(coverage * N); unique_frac ≥ stop_thresh / N
    assert (w.unique_visits >= stop_thresh).all()


# ---------------------------------------------------------------------------
# Engine-agnosticism (NX vs GT): the kernel must accept both signatures.
# ---------------------------------------------------------------------------

def test_engine_agnosticism(tmp_data_dir):
    lat_nx = Lattice2D(12, geo='sqr', pbc=True, pflip=0.1, seed=7,
                       path_data=tmp_data_dir, engine='nx')
    lat_nx.flip_random_fract_edges()
    # GT factory does not accept path_data — use default path for GT.
    lat_gt = Lattice2D(12, geo='sqr', pbc=True, pflip=0.1, seed=7,
                       engine='gt')
    # Both must expose ``.adj`` consumed by signed_lattice_tables.
    n_nx, neg_nx, fr_nx = signed_lattice_tables(lat_nx)
    n_gt, neg_gt, fr_gt = signed_lattice_tables(lat_gt)
    # Shapes must match; content may differ due to independent RNG paths.
    assert n_nx.shape == n_gt.shape == (144, 4)
    assert neg_nx.shape == neg_gt.shape
    assert fr_nx.shape == fr_gt.shape
    # Walker-on-GT also runs without error:
    rng = np.random.default_rng(0)
    start = rng.integers(0, 144, size=20, dtype=np.int64)
    out = run_walker(
        lat_gt, 'absorb', n_walkers=20, start_positions=start,
        seed=5, coverage_stop=0.2,
    )
    assert out['unique_visits'].shape == (20,)


# ---------------------------------------------------------------------------
# Kernel edge cases
# ---------------------------------------------------------------------------

def test_run_walker_rejects_bad_start_shape(disordered_lattice):
    with pytest.raises(ValueError, match="start_positions shape"):
        run_walker(disordered_lattice, 'absorb', n_walkers=10,
                   start_positions=np.array([0, 1, 2]),
                   seed=1, coverage_stop=0.2)


def test_run_walker_rejects_out_of_range_start(disordered_lattice):
    bad = np.array([9999] * 5, dtype=np.int64)
    with pytest.raises(ValueError, match="out-of-range"):
        run_walker(disordered_lattice, 'absorb', n_walkers=5,
                   start_positions=bad, seed=1, coverage_stop=0.2)


def test_run_walker_rejects_bad_rule(disordered_lattice):
    with pytest.raises(ValueError, match="rule must"):
        run_walker(disordered_lattice, 'teleport', n_walkers=5,
                   start_positions=np.zeros(5, dtype=np.int64),
                   seed=1, coverage_stop=0.2)
