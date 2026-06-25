"""Vectorized (np / cu) voter rule coverage (Phase 3).

The synchronous vectorized backend now covers the full rule family
(majority / qvoter / nonlinear), not just linear. Tests: each rule runs and
produces valid output; ``nonlinear(alpha=1)`` reproduces the linear voter; the
np backend matches the py reference statistically; the guards still reject
savedyn / the intrinsically-sequential schedules.
"""
import numpy as np
import pytest

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.VoterModel.VoterModel import VoterModel
from lrgsglib.statsys.VoterModel._vectorized_voter import CUPY_AVAILABLE


def _final_abs_magn(runlang, rule, seeds, *, side=12, steps=40, **kw):
    """Mean |final magnetization| per seed for (runlang, rule), synchronous."""
    vals = []
    for sd in seeds:
        lat = Lattice2D(side, engine="nx")
        vm = VoterModel(sg=lat, runlang=runlang, rule=rule,
                        upd_mode="synchronous", ic="random", steps=steps,
                        seed=sd, savemagn=True, savedisk=False, **kw)
        vm.run()
        vals.append(abs(float(np.mean(vm.s))))
    return np.asarray(vals)


@pytest.mark.code
@pytest.mark.parametrize("rule", ["majority", "qvoter", "nonlinear"])
def test_vectorized_rule_runs(rule):
    lat = Lattice2D(12, engine="nx")
    vm = VoterModel(sg=lat, runlang="np", rule=rule, upd_mode="synchronous",
                    ic="random", steps=30, seed=0, savemagn=True, savedisk=False)
    vm.run()
    assert set(int(x) for x in np.unique(vm.s)).issubset({-1, 1})
    assert len(vm.magn) > 0


@pytest.mark.code
def test_vectorized_qvoter_q1_is_linear_smoke():
    # q=1 q-voter copies its single sample -> the linear voter; just runs cleanly.
    lat = Lattice2D(10, engine="nx")
    vm = VoterModel(sg=lat, runlang="np", rule="qvoter", q=1, eps=0.0,
                    upd_mode="synchronous", ic="random", steps=20, seed=3,
                    savemagn=True, savedisk=False)
    vm.run()
    assert set(int(x) for x in np.unique(vm.s)).issubset({-1, 1})


@pytest.mark.physical
def test_vectorized_nonlinear_alpha1_matches_linear():
    # alpha=1 gives P(+1)=f_+ , exactly the linear voter -> same distribution.
    seeds = range(24)
    lin = _final_abs_magn("np", "linear", seeds)
    nl1 = _final_abs_magn("np", "nonlinear", seeds, alpha=1.0)
    assert abs(lin.mean() - nl1.mean()) < 0.2


@pytest.mark.physical
def test_vectorized_np_matches_py_majority():
    # py and np run the *same* synchronous majority dynamics -> close sample means.
    seeds = range(24)
    py = _final_abs_magn("py", "majority", seeds)
    npv = _final_abs_magn("np", "majority", seeds)
    assert abs(py.mean() - npv.mean()) < 0.25


@pytest.mark.physical
def test_vectorized_nonlinear_high_alpha_orders_more():
    # alpha >> 1 reinforces the local majority -> more order than the linear voter.
    seeds = range(24)
    lin = _final_abs_magn("np", "linear", seeds, steps=30)
    nl = _final_abs_magn("np", "nonlinear", seeds, steps=30, alpha=5.0)
    assert nl.mean() >= lin.mean() - 0.05


@pytest.mark.code
def test_vectorized_savedyn_still_raises():
    lat = Lattice2D(8, engine="nx")
    vm = VoterModel(sg=lat, runlang="np", rule="majority", upd_mode="synchronous",
                    steps=10, savedyn=True, savedisk=False)
    with pytest.raises(NotImplementedError):
        vm.run()


@pytest.mark.code
def test_vectorized_gillespie_still_raises():
    lat = Lattice2D(8, engine="nx")
    vm = VoterModel(sg=lat, runlang="np", rule="linear", upd_mode="gillespie",
                    steps=10, savedisk=False)
    with pytest.raises(NotImplementedError):
        vm.run()


@pytest.mark.code
@pytest.mark.skipif(not CUPY_AVAILABLE, reason="cupy not available")
def test_vectorized_cu_rule_runs():
    lat = Lattice2D(12, engine="nx")
    vm = VoterModel(sg=lat, runlang="cu", rule="nonlinear", alpha=2.0,
                    upd_mode="synchronous", ic="random", steps=20, seed=1,
                    savemagn=True, savedisk=False)
    vm.run()
    assert set(int(x) for x in np.unique(vm.s)).issubset({-1, 1})
