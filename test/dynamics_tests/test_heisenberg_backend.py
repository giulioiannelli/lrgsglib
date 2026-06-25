"""Phase-4 native parity: HeisenbergModel pybind backend + registry dispatch.

The pybind ``_heisenberg_native`` kernel wraps the SAME C
``heisenberg_metropolis_sweep`` the C-subprocess backend uses (random-node
sweeps with random-rotation proposals), so it is validated against the C
semantics and physics -- NOT bit-for-bit against the pure-Python loop, which
uses a different sampler (a random permutation per sweep) and the numpy RNG
stream.

State is an (N, 3) array of unit vectors; the kernel marshals it as a flat 3N
double buffer. Parity bar (stochastic Metropolis family):
  * valid state range (unit vectors, shape (N, 3))
  * pb determinism (same seed -> identical trajectory)
  * energy self-consistency (kernel energy == independent Python recompute)
  * physics (ordering at low T vs disorder at high T)
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys.HeisenbergModel import HeisenbergModel
from lrgsglib.statsys.HeisenbergModel.defaults import HEISENBERG_SOLVER_NAME
from lrgsglib.statsys._solver import SolverBackend
from lrgsglib.statsys._solver_engine import (
    is_backend_available,
    list_solvers_for_model,
)

pytestmark = pytest.mark.code

_PB = is_backend_available(HEISENBERG_SOLVER_NAME, SolverBackend.PB)
_needs_pb = pytest.mark.skipif(not _PB, reason="_heisenberg_native not built")


def _lattice(side=12):
    return Lattice2D(side1=side, side2=side, geo="squared", pbc=True,
                     init_nw_dict=False)


def _heisenberg(lat, *, runlang, T, steps=120, delta=0.3, seed=1):
    # NB: seed must be TRUTHY -- DynSys.__init_randomness__ does ``seed or
    # <time>``, so seed=0 is silently nondeterministic.
    m = HeisenbergModel(lat, T=T, delta=delta, steps=steps,
                        save_observables=True, runlang=runlang, seed=seed)
    m.init_heisenberg_dynamics()
    return m


def test_registry_lists_heisenberg_backends():
    backends = {b.value for b in list_solvers_for_model(HEISENBERG_SOLVER_NAME)}
    assert {"py", "pb", "c"} <= backends


@_needs_pb
def test_pb_runs_via_registry_and_state_is_valid():
    lat = _lattice()
    m = _heisenberg(lat, runlang="pb", T=0.5)
    m.run()
    assert m.s.shape == (m.N, 3)
    assert m.s.dtype == np.float64
    # All spins are unit vectors.
    norms = np.linalg.norm(m.s, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-9)
    assert len(m.ene) == m.steps and len(m.magn) == m.steps
    # Magnetisation magnitude is a valid order parameter in [0, 1].
    assert all(0.0 <= v <= 1.0 + 1e-9 for v in m.magn)


@_needs_pb
def test_pb_is_seed_reproducible():
    lat = _lattice()
    a = _heisenberg(lat, runlang="pb", T=0.5, seed=11); a.run()
    b = _heisenberg(lat, runlang="pb", T=0.5, seed=11); b.run()
    assert np.array_equal(a.s, b.s)
    assert np.allclose(a.ene, b.ene) and np.allclose(a.magn, b.magn)
    c = _heisenberg(lat, runlang="pb", T=0.5, seed=999); c.run()
    assert not np.array_equal(a.s, c.s)  # different seed -> different trajectory


@_needs_pb
def test_pb_energy_self_consistent_with_python_recompute():
    # The kernel's calc_heisenberg_energy must agree with the model's
    # independent Python energy() on the returned configuration.
    lat = _lattice()
    m = _heisenberg(lat, runlang="pb", T=0.8, seed=3)
    m.run()
    assert np.isclose(m.ene[-1], m.energy(m.s))


@_needs_pb
def test_pb_and_py_agree_on_order_disorder():
    # Different samplers, so compare the qualitative phase, not exact values.
    # On a finite 2D lattice the Heisenberg model develops strong quasi-order
    # at low T (large mean magnetisation) and stays disordered at high T.
    lat = _lattice()

    def mean_tail_order(runlang, T):
        vals = []
        for sd in (1, 2, 3):
            m = _heisenberg(lat, runlang=runlang, T=T, steps=200, seed=sd)
            m.run()
            vals.append(float(np.mean(m.magn[-40:])))
        return float(np.mean(vals))

    assert mean_tail_order("py", 0.1) > 0.5
    assert mean_tail_order("pb", 0.1) > 0.5
    assert mean_tail_order("py", 5.0) < 0.3
    assert mean_tail_order("pb", 5.0) < 0.3
