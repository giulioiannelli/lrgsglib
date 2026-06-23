"""Tests for helpers crystallized out of the VoterModel intro notebook.

``subsample`` (``utils.basic.iterables``) and the configuration-domain
observables (``edge_sign_arrays``, ``cluster_components``, ``largest_fraction``
in ``utils.statsys``) were lifted together from
``ipy/nbs/VM/VoterModel_signed_intro.ipynb``. The cluster-size *distribution* and
its native-kernel parity are exercised in depth by
``dynamics_tests/test_voter_clusters.py``; here we pin the small,
hand-verifiable behaviour of the lifted helpers themselves.
"""

import numpy as np
import pytest

from lrgsglib.utils.basic import subsample
from lrgsglib.utils.statsys import (
    cluster_components,
    edge_sign_arrays,
    largest_fraction,
)


# --------------------------------------------------------------------------- #
# subsample                                                                    #
# --------------------------------------------------------------------------- #
def test_subsample_count_and_endpoints():
    out = subsample(list(range(100)), 5)
    assert out == [0, 24, 49, 74, 99]          # endpoints kept, evenly spaced


def test_subsample_passthrough_when_short():
    assert subsample([3, 1, 4], 10) == [3, 1, 4]
    assert subsample([3, 1, 4], 3) == [3, 1, 4]   # n == k -> all of it


def test_subsample_accepts_ndarray():
    out = subsample(np.arange(10), 4)
    assert [int(x) for x in out] == [0, 3, 6, 9]


# --------------------------------------------------------------------------- #
# domain observables on a tiny 4-node path  0 - 1 - 2 - 3                      #
# --------------------------------------------------------------------------- #
def _path4_rawspin():
    idx = [np.array([1]), np.array([0, 2]), np.array([1, 3]), np.array([2])]
    b = edge_sign_arrays([np.ones(len(a), np.int8) for a in idx], "rawspin")
    return idx, b


def test_edge_sign_arrays_modes():
    signs = [np.array([-1], np.int8), np.array([-1, 1], np.int8)]
    sat = edge_sign_arrays(signs, "satisfied")
    raw = edge_sign_arrays(signs, "rawspin")
    assert [list(a) for a in sat] == [[-1], [-1, 1]]   # keeps sign(w)
    assert [list(a) for a in raw] == [[1], [1, 1]]     # all +1
    with pytest.raises(ValueError):
        edge_sign_arrays(signs, "nope")


def test_cluster_components_aligned_vs_split():
    idx, b = _path4_rawspin()
    one = cluster_components(np.array([1, 1, 1, 1], np.int8), idx, b)
    assert len(set(int(x) for x in one)) == 1          # single domain
    two = cluster_components(np.array([1, 1, -1, -1], np.int8), idx, b)
    assert len(set(int(x) for x in two)) == 2          # split at the 1-2 edge


def test_largest_fraction_on_known_configs():
    idx, b = _path4_rawspin()
    s_all = np.array([1, 1, 1, 1], np.int8)            # one domain of 4 -> 1.00
    s_half = np.array([1, 1, -1, -1], np.int8)         # two domains of 2 -> 0.50
    s_alt = np.array([1, -1, 1, -1], np.int8)          # four singletons  -> 0.25
    frac = largest_fraction([s_all, s_half, s_alt], idx, b)
    assert frac == pytest.approx([1.0, 0.5, 0.25])
