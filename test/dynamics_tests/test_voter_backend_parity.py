"""Phase-1 hardening tests for VoterModel: backend parity + API honesty.

Covers:
  * the degree-0 SIGFPE regression in the C kernel,
  * the ``save_magnetization`` contract being honoured identically by py and C,
  * the unsigned linear-voter magnetization martingale + Python<->C agreement,
  * ``upd_mode`` validation (no silently-ignored update schedules).

Cross-backend comparison is statistical, not byte-equal: the Python sampler
sweeps a random permutation while the C sampler draws N nodes with replacement,
and the C backend self-seeds from wall-clock/PID (its RNG is not yet threaded
to the Python ``seed=``).
"""

import pytest
import numpy as np


def _voter_c_binary_missing() -> bool:
    from lrgsglib.statsys.VoterModel import VoterModel
    return not (VoterModel._c_bin_dir / "VoterSimulator").exists()


# ===================================================================
# Regression: isolated (degree-0) node must not SIGFPE the C kernel
# ===================================================================


@pytest.mark.physical
@pytest.mark.integration
def test_voter_c_isolated_node_no_sigfpe(tmp_path):
    """C kernel guards ``RNG_u64() % degree`` against degree-0 nodes."""
    if _voter_c_binary_missing():
        pytest.skip("VoterSimulator binary not built")
    from lrgsglib.graphs.nx import ErdosRenyiNX
    from lrgsglib.statsys import VoterModel

    g = ErdosRenyiNX(
        n=8, p=0.0, seed=1, pflip=0.0,
        path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
    )
    assert min(dict(g.gr["G"].degree()).values()) == 0  # degree-0 node present

    v = VoterModel(sg=g, steps=10, runlang="C0", seed=1)
    v.init_voter_dynamics()
    v.run(tqdm_on=False)

    # Pre-fix this SIGFPE'd -> empty stdout -> empty self.s.
    assert v.s.size == g.N
    assert np.all(np.isin(v.s, (-1, 1)))


# ===================================================================
# save_magnetization is honoured identically by both backends
# ===================================================================


def _run_collect_magn(VoterModel, sg, steps, runlang, save):
    v = VoterModel(
        sg=sg, steps=steps, runlang=runlang, seed=1, save_magnetization=save,
    )
    v.init_voter_dynamics()
    v.run(tqdm_on=False)
    return v.magn


@pytest.mark.physical
def test_voter_save_magnetization_contract_py(tmp_path):
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    steps = 20

    def make():
        return Lattice2DNX(
            side1=6, geo="sqr", pbc=True, pflip=0.0, seed=1,
            path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
        )

    assert len(_run_collect_magn(VoterModel, make(), steps, "py", False)) == 0
    assert len(_run_collect_magn(VoterModel, make(), steps, "py", True)) == steps


@pytest.mark.physical
@pytest.mark.integration
def test_voter_save_magnetization_contract_c(tmp_path):
    if _voter_c_binary_missing():
        pytest.skip("VoterSimulator binary not built")
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    steps = 20

    def make():
        return Lattice2DNX(
            side1=6, geo="sqr", pbc=True, pflip=0.0, seed=1,
            path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
        )

    # Previously the C backend always populated magn, ignoring the flag.
    assert len(_run_collect_magn(VoterModel, make(), steps, "C0", False)) == 0
    assert len(_run_collect_magn(VoterModel, make(), steps, "C0", True)) == steps


# ===================================================================
# Unsigned linear voter: magnetization martingale + py/C agreement
# ===================================================================


@pytest.mark.physical
@pytest.mark.integration
def test_voter_martingale_and_py_c_agreement(tmp_path):
    """On a regular unsigned graph E[M_final] = M_init for both backends."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel

    N, n_up = 10, 6
    ic = np.full(N, -1, dtype=np.int8)
    ic[:n_up] = 1
    M_init = float(ic.mean())  # = 0.2
    # Statistical test: rejects the aliasing/sink bug (which drove <M> to -1)
    # while staying robust to RNG (C self-seeds from wall-clock/PID).
    R, steps, tol = 150, 100, 0.30

    g = FullyConnectedNX(
        N=N, pflip=0.0, seed=0,
        path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
    )

    def mean_final_M(runlang):
        acc = 0.0
        for r in range(R):
            v = VoterModel(sg=g, ic="custom", steps=steps, runlang=runlang, seed=r)
            v.init_voter_dynamics(custom=ic)
            v.run(tqdm_on=False)
            acc += float(np.mean(v.s))
        return acc / R

    m_py = mean_final_M("py")
    assert abs(m_py - M_init) < tol, (
        f"Python voter violates martingale: <M_final>={m_py:.3f}, M_init={M_init}"
    )

    if _voter_c_binary_missing():
        pytest.skip("VoterSimulator binary not built; Python martingale verified")
    # Both backends agreeing with theory (E[M]=M_init) implies they agree with
    # each other, without the flakier direct py-vs-C mean comparison.
    m_c = mean_final_M("C0")
    assert abs(m_c - M_init) < tol, (
        f"C voter violates martingale: <M_final>={m_c:.3f}, M_init={M_init}"
    )


# ===================================================================
# upd_mode validation (no silently-ignored schedules)
# ===================================================================


def test_voter_upd_mode_validation(tmp_path):
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0, seed=1,
        path_data=tmp_path, path_plot=tmp_path, init_nw_dict=False,
    )

    with pytest.raises(NotImplementedError):
        VoterModel(sg=lat, steps=5, runlang="py", upd_mode="synchronous")
    with pytest.raises(ValueError):
        VoterModel(sg=lat, steps=5, runlang="py", upd_mode="bogus")

    v = VoterModel(sg=lat, steps=5, runlang="py", upd_mode="asynchronous")
    assert v.upd_mode == "asynchronous"
