"""Phase-4 native parity: PottsModel pybind backend + solver-registry dispatch.

The pybind ``_potts_native`` kernel wraps the SAME C ``potts_metropolis_sweep``
the C-subprocess backend uses (random-node-with-replacement sweeps), so it is
validated against the C semantics and physics -- NOT bit-for-bit against the
pure-Python loop, which uses a different sampler (a random permutation per
sweep) and the numpy RNG stream.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys.PottsModel import PottsModel
from lrgsglib.statsys.PottsModel.defaults import POTTS_SOLVER_NAME
from lrgsglib.statsys._solver import SolverBackend
from lrgsglib.statsys._solver_engine import (
    is_backend_available,
    list_solvers_for_model,
)

pytestmark = pytest.mark.code

_PB = is_backend_available(POTTS_SOLVER_NAME, SolverBackend.PB)
_needs_pb = pytest.mark.skipif(not _PB, reason="_potts_native not built")

Q = 3


def _lattice(side=12):
    return Lattice2D(side1=side, side2=side, geo="squared", pbc=True,
                     init_nw_dict=False)


def _potts(lat, *, runlang, T, steps=120, seed=7):
    m = PottsModel(lat, q=Q, T=T, steps=steps, save_observables=True,
                   runlang=runlang, seed=seed)
    m.init_potts_dynamics()
    return m


def test_registry_lists_potts_backends():
    backends = {b.value for b in list_solvers_for_model(POTTS_SOLVER_NAME)}
    assert {"py", "pb", "c"} <= backends


@_needs_pb
def test_pb_runs_via_registry_and_state_is_valid():
    lat = _lattice()
    m = _potts(lat, runlang="pb", T=0.6)
    m.run()
    assert m.s.dtype == np.int32
    assert ((m.s >= 0) & (m.s < Q)).all()
    assert len(m.ene) == m.steps and len(m.magn) == m.steps


@_needs_pb
def test_pb_is_seed_reproducible():
    lat = _lattice()
    a = _potts(lat, runlang="pb", T=0.6, seed=11); a.run()
    b = _potts(lat, runlang="pb", T=0.6, seed=11); b.run()
    assert np.array_equal(a.s, b.s)
    assert np.allclose(a.ene, b.ene) and np.allclose(a.magn, b.magn)
    c = _potts(lat, runlang="pb", T=0.6, seed=999); c.run()
    assert not np.array_equal(a.s, c.s)  # different seed -> different trajectory


@_needs_pb
def test_pb_energy_self_consistent_with_python_recompute():
    # The kernel's calc_potts_energy must agree with the model's independent
    # Python energy() on the returned configuration.
    lat = _lattice()
    m = _potts(lat, runlang="pb", T=1.0, seed=3)
    m.run()
    assert np.isclose(m.ene[-1], m.energy(m.s))


@_needs_pb
def test_pb_and_py_agree_on_order_disorder():
    # Different samplers, so compare the qualitative phase, not exact values:
    # both order at low T and disorder at high T.
    lat = _lattice()

    def mean_tail_order(runlang, T):
        vals = []
        for sd in (1, 2, 3):
            m = _potts(lat, runlang=runlang, T=T, steps=150, seed=sd)
            m.run()
            vals.append(float(np.mean(m.magn[-30:])))
        return float(np.mean(vals))

    assert mean_tail_order("py", 0.2) > 0.5
    assert mean_tail_order("pb", 0.2) > 0.5
    assert mean_tail_order("py", 5.0) < 0.3
    assert mean_tail_order("pb", 5.0) < 0.3
