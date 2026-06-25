"""C-subprocess voter parity (Phase 2): cldist file format, seed, reproducibility.

The C subprocess now (a) seeds deterministically, (b) emits the cluster-size
distribution as a sparse-binary ``cldist`` file, decoded back to the same
``list[{size: count}]`` the py / pybind backends produce. Tests cover the file
round-trip and the C backend's new reproducibility + cluster output.
"""
import numpy as np
import pytest

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.VoterModel.VoterModel import VoterModel
from lrgsglib.statsys.VoterModel import _voter_io


# --------------------------------------------------------------------------- #
# Sparse-binary cldist file round-trip (pure, no C run)
# --------------------------------------------------------------------------- #
@pytest.mark.code
def test_load_cldist_bin_roundtrip(tmp_path):
    # two records: {size2:1, size3:1} then {size1:4}
    nrec = 2
    offsets = [0, 2, 3]          # rec0 -> entries [0,2), rec1 -> [2,3)
    sizes = [2, 3, 1]
    counts = [1, 1, 4]
    buf = np.array([nrec] + offsets + sizes + counts, dtype="<u8")
    path = tmp_path / "cldist_p=0_x.bin"
    buf.tofile(path)
    out = _voter_io.load_cldist_bin(path)
    assert out == [{2: 1, 3: 1}, {1: 4}]


@pytest.mark.code
def test_load_cldist_bin_empty(tmp_path):
    path = tmp_path / "empty.bin"
    np.array([], dtype="<u8").tofile(path)
    assert _voter_io.load_cldist_bin(path) == []


# --------------------------------------------------------------------------- #
# C subprocess: cluster emission + deterministic seed
# --------------------------------------------------------------------------- #
def _run_c(seed, *, track=True, mode="rawspin", steps=25, side=8):
    lat = Lattice2D(side, engine="nx")
    vm = VoterModel(sg=lat, runlang="C0", steps=steps, seed=seed,
                    savemagn=True, savedisk=False,
                    track_clusters=track, cluster_mode=mode)
    vm.run()
    return vm


@pytest.mark.integration
@pytest.mark.parametrize("mode", ["rawspin", "satisfied"])
def test_c_cldist_valid(mode):
    vm = _run_c(7, mode=mode)
    cd = list(vm.cluster_dist)
    assert len(cd) == 25                       # one record per recorded sweep
    # every record is a partition of the N nodes
    for rec in cd:
        assert sum(int(k) * int(v) for k, v in rec.items()) == vm.N


@pytest.mark.integration
def test_c_seed_reproducible_with_clusters():
    a, b = _run_c(11), _run_c(11)
    assert np.array_equal(a.s, b.s)
    assert np.allclose(a.magn, b.magn)
    assert list(a.cluster_dist) == list(b.cluster_dist)
    c = _run_c(12)
    assert not np.array_equal(a.s, c.s)        # a different seed diverges


@pytest.mark.integration
def test_c_seed_reproducible_magn_only():
    # reproducibility holds for plain runs too (no cluster tracking)
    a = _run_c(5, track=False)
    b = _run_c(5, track=False)
    assert np.array_equal(a.s, b.s) and np.allclose(a.magn, b.magn)


@pytest.mark.integration
def test_c_cldist_matches_py_structure_on_consensus():
    # On an unsigned lattice driven to consensus, the final record is one giant
    # cluster of all N nodes -- identical structure across backends.
    lat = Lattice2D(6, engine="nx")
    vm = VoterModel(sg=lat, runlang="C0", steps=20000, seed=1,
                    absorbing_check=True, savemagn=True, savedisk=False,
                    track_clusters=True, cluster_mode="rawspin")
    vm.run()
    assert vm.absorbed_at is not None          # reached consensus
    assert list(vm.cluster_dist)[-1] == {vm.N: 1}
