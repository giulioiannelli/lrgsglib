"""Phase-4 native parity: CoupledODEModel pybind backend + registry dispatch.

CoupledODEModel dynamics are deterministic (fixed initial condition + RHS, RK4),
so the pybind ``_code_native`` kernel and the pure-Python integrator agree
*numerically* (not merely statistically) -- the strongest parity bar in the
Phase-4 model set. Mirrors ``test_kuramoto_backend.py``.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys.CoupledODEModel import CoupledODEModel
from lrgsglib.statsys.CoupledODEModel.defaults import CODE_SOLVER_NAME
from lrgsglib.statsys._solver import SolverBackend
from lrgsglib.statsys._solver_engine import (
    is_backend_available,
    list_solvers_for_model,
)

pytestmark = pytest.mark.code

_PB = is_backend_available(CODE_SOLVER_NAME, SolverBackend.PB)
_needs_pb = pytest.mark.skipif(not _PB, reason="_code_native not built")


def _lattice(side=10, pflip=0.0):
    return Lattice2D(side1=side, side2=side, geo="squared", pbc=True,
                     pflip=pflip, init_nw_dict=False)


def _code(lat, *, runlang, coupling_type="linear", local_type="fitzhugh",
          coupling_strength=0.5, steps=200, seed=1):
    m = CoupledODEModel(lat, coupling_type=coupling_type, local_type=local_type,
                        coupling_strength=coupling_strength, dt=0.01,
                        steps=steps, runlang=runlang, seed=seed)
    m.init_ode_dynamics()
    return m


def test_registry_lists_code_backends():
    backends = {b.value for b in list_solvers_for_model(CODE_SOLVER_NAME)}
    assert {"py", "pb", "c"} <= backends


@_needs_pb
def test_pb_runs_via_registry():
    lat = _lattice()
    m = _code(lat, runlang="pb")
    m.run()
    assert m.s.dtype == np.float64 and m.s.shape == (lat.N,)
    assert np.all(np.isfinite(m.s))


@_needs_pb
def test_pb_matches_python_numerically_linear():
    # Deterministic dynamics on an UNSIGNED lattice (pflip=0): linear/diffusive
    # coupling sums w*(x_j - x_i); with all-positive weights the C kernel's
    # signed degree sum(w) equals the Python abs degree sum(|w|), so py and pb
    # trajectories coincide to machine precision.
    lat = _lattice(pflip=0.0)
    mpy = _code(lat, runlang="py", coupling_type="linear"); mpy.run()
    mpb = _code(lat, runlang="pb", coupling_type="linear"); mpb.run()
    assert np.allclose(mpy.sini, mpb.sini)  # same initial condition
    d = np.max(np.abs(np.asarray(mpy.s) - np.asarray(mpb.s)))
    assert d < 1e-9, f"linear-coupling state mismatch {d:.3e}"


@_needs_pb
def test_pb_product_signed_finite_and_deterministic():
    # The native kernel mirrors the C-subprocess product coupling g = w*x_j*x_i.
    # NOTE: the model's *vectorised* Python path (_g_coupling_vec) implements a
    # different product convention (sum_j A_ij x_j^2 vs the C kernel's
    # sum_j A_ij x_j x_i), a pre-existing py-vs-C divergence -- so we do not
    # assert py<->pb numerical equality for product coupling. Instead we check
    # the pb backend stays finite on a SIGNED (frustrated) lattice and is
    # bit-reproducible for a fixed IC.
    lat = _lattice(pflip=0.3)
    rng = np.random.default_rng(1)
    x0 = rng.uniform(-0.1, 0.1, size=lat.N).astype(np.float64)
    common = dict(coupling_type="product", local_type="none",
                  coupling_strength=0.05, dt=0.01, steps=50, ic="custom", seed=1)
    m1 = CoupledODEModel(lat, runlang="pb", **common)
    m1.init_ode_dynamics(custom=x0); m1.run()
    m2 = CoupledODEModel(lat, runlang="pb", **common)
    m2.init_ode_dynamics(custom=x0); m2.run()
    assert np.all(np.isfinite(m1.s))
    assert np.array_equal(m1.s, m2.s)


@_needs_pb
def test_pb_deterministic_same_seed():
    # Same truthy seed -> identical initial condition -> identical pb trajectory.
    lat = _lattice()
    m1 = _code(lat, runlang="pb", seed=1); m1.run()
    m2 = _code(lat, runlang="pb", seed=1); m2.run()
    assert np.array_equal(m1.s, m2.s)


@_needs_pb
def test_pb_matches_python_logistic():
    # Local logistic term with explicit (default-equal) params; still numerically
    # exact on an unsigned lattice with linear coupling.
    lat = _lattice(pflip=0.0)
    kw = dict(coupling_type="linear", local_type="logistic",
              coupling_strength=0.2, steps=150)
    mpy = _code(lat, runlang="py", **kw); mpy.run()
    mpb = _code(lat, runlang="pb", **kw); mpb.run()
    d = np.max(np.abs(np.asarray(mpy.s) - np.asarray(mpb.s)))
    assert d < 1e-9, f"logistic state mismatch {d:.3e}"
