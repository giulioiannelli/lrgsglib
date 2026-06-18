"""Phase-2 tests for the in-process pybind11 voter backend (``pb_voter``).

The pybind backend reuses the same ``voter_model_Nstep`` C kernel as the
subprocess backend but passes the graph as numpy CSR arrays: no file I/O, GT
graphs supported, and -- unlike the subprocess backend -- seed-reproducible.
"""

import pytest
import numpy as np


def _pb_available() -> bool:
    try:
        from lrgsglib.statsys.VoterModel.ccore import _voter_native  # noqa: F401
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _pb_available(), reason="_voter_native pybind module not built"
)


@pytest.mark.physical
def test_pb_voter_runs_and_magn_contract(tmp_path):
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    def make():
        return Lattice2DNX(
            side1=12, geo="sqr", pbc=True, pflip=0.0, seed=7,
            path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
        )

    v = VoterModel(sg=make(), steps=40, runlang="pb_voter", seed=7,
                   save_magnetization=True)
    v.run(tqdm_on=False)
    assert v.s.size == v.N
    assert np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) == 40

    # save_magnetization=False must leave magn empty (parity with py/C)
    v2 = VoterModel(sg=make(), steps=40, runlang="pb_voter", seed=7,
                    save_magnetization=False)
    v2.run(tqdm_on=False)
    assert len(v2.magn) == 0


def test_pb_voter_seed_reproducible(tmp_path):
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    def run(seed):
        lat = Lattice2DNX(
            side1=10, geo="sqr", pbc=True, pflip=0.0, seed=seed,
            path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
        )
        m = VoterModel(sg=lat, steps=60, runlang="pb_voter", seed=seed)
        m.run(tqdm_on=False)
        return m.s.copy()

    assert np.array_equal(run(2024), run(2024))  # same seed -> identical run


@pytest.mark.physical
def test_pb_voter_martingale(tmp_path):
    """Unsigned linear voter conserves E[M] under the pybind backend too."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel

    N, n_up = 10, 6
    ic = np.full(N, -1, dtype=np.int8)
    ic[:n_up] = 1
    M_init = float(ic.mean())  # 0.2
    R, tol = 200, 0.30

    g = FullyConnectedNX(
        N=N, pflip=0.0, seed=1,
        path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
    )
    acc = 0.0
    for r in range(1, R + 1):
        m = VoterModel(sg=g, ic="custom", steps=100, runlang="pb_voter", seed=r)
        m.init_voter_dynamics(custom=ic)
        m.run(tqdm_on=False)
        acc += float(np.mean(m.s))
    assert abs(acc / R - M_init) < tol, f"<M_final>={acc / R:.3f}, M_init={M_init}"


@pytest.mark.integration
def test_pb_voter_on_gt_graph():
    """The pybind backend runs on a graph-tool graph (no file I/O)."""
    pytest.importorskip("graph_tool")
    from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT
    from lrgsglib.statsys import VoterModel

    g = Lattice2DGT(side1=12, geo="sqr", periodic=True, pflip=0.0, seed=3)
    v = VoterModel(sg=g, steps=50, runlang="pb_voter", seed=3,
                   save_magnetization=True)
    v.run(tqdm_on=False)
    assert v.s.size == g.N
    assert np.all(np.isin(v.s, (-1, 1)))
    assert len(v.magn) == 50
