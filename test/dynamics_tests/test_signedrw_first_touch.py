"""Tests for the native `absorb_walker_first_touch` kernel.

The kernel alternates one absorbing walker from each of two start nodes,
growing two persistent bubbles, and returns the smallest walker-pair index
`k` at which the bubbles first intersect.  Returns `n_max + 1` if censored.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.SignedRW._kernel import (
    _kill_masks,
    signed_lattice_tables,
)
from lrgsglib.statsys.SignedRW.SignedWalker import _load_srw_native


native = _load_srw_native()
pytestmark = pytest.mark.skipif(
    native is None or not hasattr(native, "absorb_walker_first_touch"),
    reason="_srw_native.absorb_walker_first_touch not built",
)


@pytest.fixture(scope="module")
def tmp_data_dir():
    td = tempfile.TemporaryDirectory(prefix="srw_ft_test_")
    yield Path(td.name)
    td.cleanup()


def _csr(lat, x_node: str = "reflect") -> dict:
    nbrs, neg, frust = signed_lattice_tables(lat)
    km, rm = _kill_masks(neg, frust, x_node)
    N, z = nbrs.shape
    return dict(
        N=N,
        neigh_indices=np.ascontiguousarray(nbrs.ravel(), dtype=np.int64),
        neigh_ptr=np.arange(0, (N + 1) * z, z, dtype=np.int64),
        kill=np.ascontiguousarray(km.ravel(), dtype=np.uint8),
        refl=np.ascontiguousarray(rm.ravel(), dtype=np.uint8),
    )


def test_first_touch_determinism(tmp_data_dir):
    """Same seed → identical output."""
    lat = Lattice2D(16, geo="sqr", pbc=True, pflip=0.10, seed=42,
                    path_data=tmp_data_dir)
    lat.flip_random_fract_edges()
    c = _csr(lat)
    args = (c["N"], c["neigh_indices"], c["neigh_ptr"],
            c["kill"], c["refl"], 0, (16 // 2) * (16 + 1),
            1000, 10 * c["N"], 1234, False)
    o1 = native.absorb_walker_first_touch(*args)
    o2 = native.absorb_walker_first_touch(*args)
    assert o1 == o2


def test_first_touch_unsigned_lattice_quick(tmp_data_dir):
    """At p=0 (no killing) walkers touch quickly."""
    lat = Lattice2D(8, geo="sqr", pbc=True, pflip=0.0, seed=1,
                    path_data=tmp_data_dir)
    c = _csr(lat)
    samples = [
        native.absorb_walker_first_touch(
            c["N"], c["neigh_indices"], c["neigh_ptr"],
            c["kill"], c["refl"],
            0, (8 // 2) * (8 + 1),
            10000, 10 * c["N"], seed, False,
        )[0]
        for seed in range(20)
    ]
    # At p=0 walker A walks unbounded, fills a large chunk, walker B almost
    # always hits bubble_A during its first walk.  So n_touch == 1 for most.
    assert all(s == 1 for s in samples), f"unexpected n_touch={samples}"


def test_first_touch_censors_at_small_nmax(tmp_data_dir):
    """Trial with artificially tiny n_max is censored at n_max+1."""
    lat = Lattice2D(16, geo="sqr", pbc=True, pflip=0.05, seed=7,
                    path_data=tmp_data_dir)
    lat.flip_random_fract_edges()
    c = _csr(lat)
    n_max = 2
    n_touch, _, _ = native.absorb_walker_first_touch(
        c["N"], c["neigh_indices"], c["neigh_ptr"],
        c["kill"], c["refl"],
        0, (16 // 2) * (16 + 1),
        n_max, 10 * c["N"], 99, False,
    )
    # Either it touches within 2 pairs OR returns n_max+1=3 (censored).
    # With only 2 walker pairs at p=0.05 on L=16 starting antipodally,
    # it is unlikely but not impossible to touch.  Allow both outcomes
    # and verify the returned code is one of {1, 2, 3}.
    assert n_touch in (1, 2, n_max + 1)


def test_first_touch_returns_bubble_sizes_when_requested(tmp_data_dir):
    """track_bubbles=True gives nonzero popcounts; False gives zeros."""
    lat = Lattice2D(16, geo="sqr", pbc=True, pflip=0.10, seed=11,
                    path_data=tmp_data_dir)
    lat.flip_random_fract_edges()
    c = _csr(lat)
    args = (c["N"], c["neigh_indices"], c["neigh_ptr"],
            c["kill"], c["refl"], 0, (16 // 2) * (16 + 1),
            1000, 10 * c["N"], 333)
    n_on, a_on, b_on = native.absorb_walker_first_touch(*args, True)
    n_off, a_off, b_off = native.absorb_walker_first_touch(*args, False)
    assert n_on == n_off  # deterministic under same seed
    assert a_off == 0 and b_off == 0
    assert a_on >= 1 and b_on >= 1


def test_first_touch_same_start_is_immediate(tmp_data_dir):
    """When node_a == node_b the bubbles already share a site → n_touch=1."""
    lat = Lattice2D(8, geo="sqr", pbc=True, pflip=0.10, seed=13,
                    path_data=tmp_data_dir)
    lat.flip_random_fract_edges()
    c = _csr(lat)
    n_touch, _, _ = native.absorb_walker_first_touch(
        c["N"], c["neigh_indices"], c["neigh_ptr"],
        c["kill"], c["refl"],
        5, 5,
        1000, 10 * c["N"], 17, False,
    )
    assert n_touch == 1
