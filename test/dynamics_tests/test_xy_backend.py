"""Phase-4 native parity: XYModel pybind backend + solver-registry dispatch.

The pybind ``_xy_native`` kernel wraps the SAME C ``xy_metropolis_sweep`` the
C-subprocess backend uses (random-node-with-replacement sweeps), so it is
validated against the C semantics and physics -- NOT bit-for-bit against the
pure-Python loop, which uses a different sampler (a random permutation per
sweep) and the numpy RNG stream.

XYModel is a Metropolis VecDynSys with a float64 angle (radians) per node, so
the parity bar is statistical (ordered vs disordered phase, valid state range)
plus pb determinism (same seed -> identical) and kernel-energy self-consistency
(the recorded energy must match an independent Python recompute).
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys.XYModel import XYModel
from lrgsglib.statsys.XYModel.defaults import XY_SOLVER_NAME
from lrgsglib.statsys._solver import SolverBackend
from lrgsglib.statsys._solver_engine import (
    is_backend_available,
    list_solvers_for_model,
)

pytestmark = pytest.mark.code

_PB = is_backend_available(XY_SOLVER_NAME, SolverBackend.PB)
_needs_pb = pytest.mark.skipif(not _PB, reason="_xy_native not built")


def _lattice(side=12):
    return Lattice2D(side1=side, side2=side, geo="squared", pbc=True,
                     init_nw_dict=False)


def _xy(lat, *, runlang, T, steps=120, seed=7):
    # NOTE: seed must be TRUTHY -- DynSys.__init_randomness__ does
    # ``seed or <time>``, so seed=0 is silently nondeterministic.
    m = XYModel(lat, T=T, delta=0.5, steps=steps, save_observables=True,
                runlang=runlang, seed=seed)
    m.init_xy_dynamics()
    return m


def test_registry_lists_xy_backends():
    backends = {b.value for b in list_solvers_for_model(XY_SOLVER_NAME)}
    assert {"py", "pb", "c"} <= backends


@_needs_pb
def test_pb_runs_via_registry_and_state_is_valid():
    lat = _lattice()
    m = _xy(lat, runlang="pb", T=0.6)
    m.run()
    assert m.s.dtype == np.float64
    assert m.s.shape == (m.N,)
    # Angles live in [0, 2*pi) (proposal is wrapped into that range).
    assert np.all(np.isfinite(m.s))
    assert m.s.min() >= 0.0 and m.s.max() < 2.0 * np.pi + 1e-9
    assert len(m.ene) == m.steps and len(m.magn) == m.steps
    # Magnetisation is |<exp(i theta)>| in [0, 1].
    assert all(0.0 <= x <= 1.0 + 1e-9 for x in m.magn)


@_needs_pb
def test_pb_is_seed_reproducible():
    lat = _lattice()
    a = _xy(lat, runlang="pb", T=0.6, seed=11); a.run()
    b = _xy(lat, runlang="pb", T=0.6, seed=11); b.run()
    assert np.array_equal(a.s, b.s)
    assert np.allclose(a.ene, b.ene) and np.allclose(a.magn, b.magn)
    c = _xy(lat, runlang="pb", T=0.6, seed=999); c.run()
    assert not np.array_equal(a.s, c.s)  # different seed -> different trajectory


@_needs_pb
def test_pb_energy_self_consistent_with_python_recompute():
    # The kernel's calc_xy_energy must agree with the model's independent
    # Python energy() on the returned configuration.
    lat = _lattice()
    m = _xy(lat, runlang="pb", T=1.0, seed=3)
    m.run()
    assert np.isclose(m.ene[-1], m.energy(m.s), rtol=1e-9, atol=1e-9)


@_needs_pb
def test_pb_and_py_agree_on_order_disorder():
    # Different samplers, so compare the qualitative phase, not exact values:
    # the XY magnetisation orders (large |m|) at low T and disorders (small |m|)
    # at high T. Average over a few seeds + the trajectory tail to damp noise.
    lat = _lattice()

    def mean_tail_order(runlang, T):
        vals = []
        for sd in (1, 2, 3):
            m = _xy(lat, runlang=runlang, T=T, steps=200, seed=sd)
            m.run()
            vals.append(float(np.mean(m.magn[-40:])))
        return float(np.mean(vals))

    assert mean_tail_order("py", 0.05) > 0.5
    assert mean_tail_order("pb", 0.05) > 0.5
    assert mean_tail_order("py", 5.0) < 0.35
    assert mean_tail_order("pb", 5.0) < 0.35
