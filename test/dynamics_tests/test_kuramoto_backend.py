"""Phase-4 native parity: KuramotoModel pybind backend + registry dispatch.

Kuramoto dynamics are deterministic (fixed omega + initial phases, RK4), so the
pybind ``_kuramoto_native`` kernel and the pure-Python integrator agree
*numerically* (not merely statistically) -- the strongest parity bar in the
Phase-4 model set.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys.KuramotoModel import KuramotoModel
from lrgsglib.statsys.KuramotoModel.defaults import KURAMOTO_SOLVER_NAME
from lrgsglib.statsys._solver import SolverBackend
from lrgsglib.statsys._solver_engine import (
    is_backend_available,
    list_solvers_for_model,
)

pytestmark = pytest.mark.code

_PB = is_backend_available(KURAMOTO_SOLVER_NAME, SolverBackend.PB)
_needs_pb = pytest.mark.skipif(not _PB, reason="_kuramoto_native not built")


def _lattice(side=10):
    return Lattice2D(side1=side, side2=side, geo="squared", pbc=True,
                     init_nw_dict=False)


def _kuramoto(lat, *, runlang, coupling=2.0, steps=200, seed=1):
    omega = np.linspace(-0.5, 0.5, lat.N)  # deterministic, explicit
    m = KuramotoModel(lat, coupling=coupling, omega=omega, dt=0.01,
                      steps=steps, save_order_param=True, runlang=runlang,
                      seed=seed)
    m.init_kuramoto_dynamics()
    return m


def test_registry_lists_kuramoto_backends():
    backends = {b.value for b in list_solvers_for_model(KURAMOTO_SOLVER_NAME)}
    assert {"py", "pb", "c"} <= backends


@_needs_pb
def test_pb_runs_via_registry():
    lat = _lattice()
    m = _kuramoto(lat, runlang="pb")
    m.run()
    assert m.s.dtype == np.float64 and m.s.shape == (lat.N,)
    assert ((m.s >= 0.0) & (m.s < 2 * np.pi + 1e-9)).all()  # wrapped to [0, 2pi)
    assert len(m.order_params) == m.steps


@_needs_pb
def test_pb_matches_python_numerically():
    # Deterministic dynamics -> py and pb trajectories must coincide to RK4
    # tolerance (here, essentially machine precision).
    lat = _lattice()
    mpy = _kuramoto(lat, runlang="py"); mpy.run()
    mpb = _kuramoto(lat, runlang="pb"); mpb.run()
    assert np.allclose(mpy.sini, mpb.sini)  # same initial condition
    d_order = np.max(np.abs(np.array(mpy.order_params) - np.array(mpb.order_params)))
    assert d_order < 1e-9, f"order-parameter mismatch {d_order:.3e}"
    dtheta = np.abs((np.asarray(mpy.s) - np.asarray(mpb.s) + np.pi) % (2 * np.pi) - np.pi)
    assert np.max(dtheta) < 1e-8, f"phase mismatch {np.max(dtheta):.3e}"


@_needs_pb
def test_pb_coupling_increases_coherence():
    # Identical natural frequencies: stronger coupling drives more coherence
    # (monotone in K). Robust to the mean-field K/N normalization on a sparse
    # lattice, where absolute synchronization needs a long integration time.
    lat = _lattice()
    omega = np.zeros(lat.N)

    def final_order(K):
        m = KuramotoModel(lat, coupling=K, omega=omega, dt=0.01, steps=600,
                          save_order_param=True, runlang="pb", seed=1)
        m.init_kuramoto_dynamics()
        m.run()
        return m.order_params[-1]

    assert final_order(50.0) > final_order(0.5)
