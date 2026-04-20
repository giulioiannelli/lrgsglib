"""Unit tests for :mod:`lrgsglib.statsys.SignedRW._overlap`."""

from __future__ import annotations

import numpy as np
import pytest

from lrgsglib.statsys.SignedRW._overlap import (
    aggregate_visits,
    overlap,
    overlap_all,
    pairwise_overlap_matrix,
)


# ---------------------------------------------------------------------------
# aggregate_visits
# ---------------------------------------------------------------------------

def test_aggregate_sums_per_walker_counts():
    vw = np.array([[1, 2, 0, 0],
                   [0, 1, 1, 0],
                   [0, 0, 3, 2]], dtype=np.int32)
    np.testing.assert_array_equal(
        aggregate_visits(vw),
        np.array([1, 3, 4, 2], dtype=np.int64),
    )


def test_aggregate_rejects_1d_input():
    with pytest.raises(ValueError, match="shape"):
        aggregate_visits(np.array([1, 2, 3]))


# ---------------------------------------------------------------------------
# overlap — core identities
# ---------------------------------------------------------------------------

def test_l2_raw_equals_dot_product():
    c1 = np.array([2, 0, 3, 1], dtype=np.int64)
    c2 = np.array([1, 4, 0, 2], dtype=np.int64)
    assert overlap(c1, c2, kind='l2_raw') == pytest.approx(2 + 0 + 0 + 2)


def test_l2_raw_self_equals_sum_of_squares():
    c = np.array([2, 0, 3, 1], dtype=np.int64)
    assert overlap(c, c, kind='l2_raw') == pytest.approx(4 + 0 + 9 + 1)


def test_l2_norm_is_bounded_by_maxnorm():
    c = np.array([1, 2, 3, 4], dtype=np.int64)
    # ⟨ρ, ρ⟩ = Σ ρ² and Σ ρ = 1, so ⟨ρ,ρ⟩ ≤ 1
    o = overlap(c, c, kind='l2_norm')
    assert 0 <= o <= 1


def test_l2_rescaled_is_N_times_l2_norm():
    c1 = np.array([1, 2, 3, 4], dtype=np.int64)
    c2 = np.array([4, 3, 2, 1], dtype=np.int64)
    N = c1.size
    a = overlap(c1, c2, kind='l2_rescaled')
    b = N * overlap(c1, c2, kind='l2_norm')
    assert a == pytest.approx(b)


def test_fidelity_self_is_one():
    c = np.array([2, 0, 3, 1], dtype=np.int64)
    assert overlap(c, c, kind='fidelity') == pytest.approx(1.0)


def test_fidelity_is_in_unit_interval():
    c1 = np.array([1, 2, 3, 4], dtype=np.int64)
    c2 = np.array([4, 3, 2, 1], dtype=np.int64)
    f = overlap(c1, c2, kind='fidelity')
    assert 0 <= f <= 1 + 1e-12


def test_tv_self_is_zero():
    c = np.array([2, 0, 3, 1], dtype=np.int64)
    assert overlap(c, c, kind='tv') == pytest.approx(0.0)


def test_tv_disjoint_supports_equals_one():
    c1 = np.array([3, 0, 0, 0], dtype=np.int64)
    c2 = np.array([0, 0, 0, 5], dtype=np.int64)
    assert overlap(c1, c2, kind='tv') == pytest.approx(1.0)


def test_jaccard_disjoint_supports_is_zero():
    c1 = np.array([3, 0, 0, 0], dtype=np.int64)
    c2 = np.array([0, 0, 0, 5], dtype=np.int64)
    assert overlap(c1, c2, kind='jaccard') == 0.0


def test_jaccard_identical_supports_is_one():
    c1 = np.array([2, 5, 1, 0], dtype=np.int64)
    c2 = np.array([7, 4, 3, 0], dtype=np.int64)
    assert overlap(c1, c2, kind='jaccard') == 1.0


def test_zero_field_yields_zero_density_overlap():
    c = np.array([2, 0, 3, 1], dtype=np.int64)
    z = np.zeros_like(c)
    # raw L² with zero is 0
    assert overlap(c, z, kind='l2_raw') == 0.0
    # density-based kinds use _as_density which maps zero-sum → zero vector
    assert overlap(c, z, kind='l2_norm') == 0.0
    assert overlap(c, z, kind='fidelity') == 0.0


# ---------------------------------------------------------------------------
# overlap — input-validation
# ---------------------------------------------------------------------------

def test_bad_kind_raises():
    c = np.array([1, 2, 3], dtype=np.int64)
    with pytest.raises(ValueError, match="kind must"):
        overlap(c, c, kind='entropy')


def test_shape_mismatch_raises():
    with pytest.raises(ValueError, match="shape mismatch"):
        overlap(np.array([1, 2, 3]), np.array([1, 2]))


def test_negative_counts_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        overlap(np.array([-1, 0, 0]), np.array([1, 0, 0]))


# ---------------------------------------------------------------------------
# overlap_all / pairwise
# ---------------------------------------------------------------------------

def test_overlap_all_covers_every_kind():
    c1 = np.array([1, 2, 3], dtype=np.int64)
    c2 = np.array([4, 2, 1], dtype=np.int64)
    d = overlap_all(c1, c2)
    assert set(d) == {'l2_raw', 'l2_norm', 'l2_rescaled',
                      'fidelity', 'tv', 'jaccard'}
    for v in d.values():
        assert np.isfinite(v)


def test_pairwise_matrix_is_symmetric():
    vw = np.array([[1, 2, 0, 0],
                   [0, 1, 1, 0],
                   [0, 0, 3, 2]], dtype=np.int64)
    m = pairwise_overlap_matrix(vw, kind='l2_raw')
    assert m.shape == (3, 3)
    np.testing.assert_allclose(m, m.T)
    # diagonal equals Σ vw[i]²
    for i in range(3):
        assert m[i, i] == pytest.approx(float((vw[i].astype(np.int64) ** 2).sum()))
