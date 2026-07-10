"""Phase-4 native parity: ReactionDiffusionModel pybind backend + registry dispatch.

Reaction-diffusion dynamics are deterministic (fixed diffusion + reaction +
initial field, RK4 with per-step clamp to [0, inf)), so the pybind ``_rd_native``
kernel and the pure-Python integrator agree *numerically* (not merely
statistically) -- the strongest parity bar in the Phase-4 model set.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys._solver import SolverBackend
from lrgsglib.statsys._solver_engine import (
    is_backend_available,
    list_solvers_for_model,
)
from lrgsglib.statsys.ReactionDiffusionModel import ReactionDiffusionModel
from lrgsglib.statsys.ReactionDiffusionModel.defaults import RD_SOLVER_NAME

pytestmark = pytest.mark.code

_PB = is_backend_available(RD_SOLVER_NAME, SolverBackend.PB)
_needs_pb = pytest.mark.skipif(not _PB, reason="_rd_native not built")


def _lattice(tmp_path, side=10):
    # path_data=tmp_path: savedisk outputs must never land in the repo data/
    return Lattice2D(
        side1=side,
        side2=side,
        geo="squared",
        pbc=True,
        init_nw_dict=False,
        path_data=tmp_path,
    )


def _rd(
    lat,
    *,
    runlang,
    reaction="fisher_kpp",
    reaction_params=None,
    D=0.1,
    dt=0.01,
    steps=200,
    seed=1,
    ic="uniform",
):
    # Deterministic, explicit initial field (no RNG dependence on seed=0 gotcha).
    u0 = np.linspace(0.1, 0.9, lat.N)
    m = ReactionDiffusionModel(
        lat,
        D=D,
        reaction=reaction,
        reaction_params=reaction_params or {"r": 1.0},
        dt=dt,
        steps=steps,
        runlang=runlang,
        seed=seed,
        ic="custom",
    )
    m.init_rd_dynamics(custom=u0)
    return m


def test_registry_lists_rd_backends():
    backends = {b.value for b in list_solvers_for_model(RD_SOLVER_NAME)}
    assert {"py", "pb", "c"} <= backends


@_needs_pb
def test_pb_runs_via_registry(tmp_path):
    lat = _lattice(tmp_path)
    m = _rd(lat, runlang="pb")
    m.run()
    assert m.s.dtype == np.float64 and m.s.shape == (lat.N,)
    assert (m.s >= 0.0).all()  # clamped to non-negative concentrations


@_needs_pb
@pytest.mark.parametrize(
    "reaction,params",
    [
        ("fisher_kpp", {"r": 1.0}),
        ("bistable", {"a": 0.3}),
        ("linear", {"r": 0.5}),
        ("none", {}),
    ],
)
def test_pb_matches_python_numerically(reaction, params, tmp_path):
    # Deterministic dynamics -> py and pb final fields must coincide to RK4
    # tolerance (essentially machine precision).
    lat = _lattice(tmp_path)
    mpy = _rd(lat, runlang="py", reaction=reaction, reaction_params=params)
    mpb = _rd(lat, runlang="pb", reaction=reaction, reaction_params=params)
    mpy.run()
    mpb.run()
    assert np.allclose(mpy.sini, mpb.sini)  # same initial condition
    d_state = np.max(np.abs(np.asarray(mpy.s) - np.asarray(mpb.s)))
    assert d_state < 1e-9, f"[{reaction}] final-field mismatch {d_state:.3e}"
    # Observable parity: mean concentration of the final field.
    d_mean = abs(mpy.mean_concentration() - mpb.mean_concentration())
    assert (
        d_mean < 1e-9
    ), f"[{reaction}] mean-concentration mismatch {d_mean:.3e}"


@_needs_pb
def test_pb_deterministic(tmp_path):
    # Same inputs -> identical output (kernel is deterministic).
    lat = _lattice(tmp_path)
    m1 = _rd(lat, runlang="pb")
    m1.run()
    m2 = _rd(lat, runlang="pb")
    m2.run()
    assert np.array_equal(m1.s, m2.s)


@_needs_pb
def test_pb_fisher_kpp_grows_toward_carrying_capacity(tmp_path):
    # Fisher-KPP logistic growth: a sub-capacity field (u < 1) must increase
    # its mean concentration over time (toward the u=1 stable state).
    lat = _lattice(tmp_path)
    m = _rd(
        lat,
        runlang="pb",
        reaction="fisher_kpp",
        reaction_params={"r": 1.0},
        D=0.05,
        steps=400,
    )
    mean0 = m.mean_concentration()
    m.run()
    assert m.mean_concentration() > mean0
    assert (m.s <= 1.0 + 1e-6).all()  # stays below carrying capacity
