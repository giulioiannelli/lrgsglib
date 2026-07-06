"""XY vertical slice: XYBase + XYMetropolis on the shared spine.

The third consumer of the statsys-unification architecture (D-B8/D-B9, XY
after Ising and Potts) and the FIRST continuous-state substrate: the SAME
``statsys._thermal`` engines drive float64 angles — no engine code changed.
Pins:

- the §3b Hamiltonian bookkeeping (ΔE↔E from ONE symmetric J; incremental
  running energy == full recompute after every run);
- the exact Ising limit (angles restricted to {0, π} give
  ``cos(θ_i − θ_j) = s_i s_j``, so XY energies and π-rotation ΔE's equal the
  Ising ones on the same graph);
- physics sanity per move/update mode (quench monotonicity, order/disorder
  around the BKT scale, angle-multiset conservation under Kawasaki, cluster
  ordering — including on SIGNED graphs, where the embedded-Ising
  reflection keeps cluster moves valid, ref M6);
- the capability boundaries (heatbath/gillespie not wired for a continuous
  state, tie_flip_p meaningless, pb kernel guards) — hard errors,
  invariant #3;
- the legacy XYModel stays untouched and agrees on the Hamiltonian.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import IsingMetropolis
from lrgsglib.statsys.XYModel import XYMetropolis, XYModel

SIDE = 8
N = SIDE * SIDE
SEED = 42


def make_lattice(tmp_path, pflip=0.1):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=7,
    )


def _native_available() -> bool:
    try:
        from lrgsglib.statsys.XYModel.ccore import _xy_native  # noqa: F401

        return True
    except Exception:
        return False


needs_native = pytest.mark.skipif(
    not _native_available(), reason="_xy_native not built"
)


# =========================================================================
# Bookkeeping: the §3b contract on the continuous substrate
# =========================================================================
def test_defaults_construct_and_run(tmp_path):
    sg = make_lattice(tmp_path)
    model = XYMetropolis(sg, T=1.0, steps=5, seed=SEED, savedisk=False)
    model.run()
    assert len(model.ene) == 5 and len(model.magn) == 5
    assert model.s.dtype == np.float64
    assert ((model.s >= 0.0) & (model.s < 2.0 * np.pi)).all()
    # record-then-sweep cadence + incremental energy invariant
    assert np.isclose(model.ene[0], model.energy_intensive(model.sini))
    assert np.isclose(model._e_running, model.total_energy())


def test_hamiltonian_matches_legacy(tmp_path):
    """The new CSR total_energy equals the legacy dense-loop energy()."""
    sg = make_lattice(tmp_path, pflip=0.2)
    new = XYMetropolis(sg, T=1.0, steps=1, seed=SEED, savedisk=False)
    new.init_dynamics()
    legacy = XYModel(sg, T=1.0, steps=1, seed=SEED)
    assert np.isclose(new.total_energy(), legacy.energy(new.s))


def test_binary_angles_are_ising(tmp_path):
    """Angles restricted to {0, π} give cos(θ_i − θ_j) = s_i s_j under
    s = cos(θ), so the XY Hamiltonian and the π-rotation ΔE equal the Ising
    ones on the same signed graph, exactly."""
    sg = make_lattice(tmp_path, pflip=0.15)
    rng = np.random.RandomState(3)
    spins = rng.choice([-1, 1], size=sg.N).astype(np.int8)
    angles = np.where(spins > 0, 0.0, np.pi)
    xy = XYMetropolis(sg, T=1.0, steps=1, seed=3, ic="custom", savedisk=False)
    xy.init_dynamics(custom=angles)
    ising = IsingMetropolis(sg, T=1.0, steps=1, seed=3, savedisk=False)
    ising.init_dynamics()
    ising.s = spins
    assert np.isclose(xy.total_energy(), ising.total_energy())
    for nd in range(sg.N):
        flipped = float((angles[nd] + np.pi) % (2.0 * np.pi))
        dE_xy = xy.delta_E(nd, flipped)
        dE_is = ising.delta_E(nd, -int(spins[nd]))
        assert np.isclose(dE_xy, dE_is)


def test_seeded_reproducibility(tmp_path):
    finals = []
    for _ in range(2):
        sg = make_lattice(tmp_path)
        m = XYMetropolis(sg, T=1.5, steps=6, seed=SEED, savedisk=False)
        m.run()
        finals.append(m.s.copy())
    assert np.array_equal(finals[0], finals[1])


# =========================================================================
# Physics: quench, order/disorder, update modes
# =========================================================================
def test_zero_T_quench_is_monotone(tmp_path):
    """At T = 0 only ΔE <= 0 proposals are accepted (ΔE = 0 has measure
    zero), so the energy trace is nonincreasing."""
    sg = make_lattice(tmp_path, pflip=0.1)
    m = XYMetropolis(sg, T=0.0, steps=20, seed=SEED, savedisk=False)
    m.run()
    assert np.all(np.diff(np.asarray(m.ene)) <= 1e-12)


def test_orders_at_low_T_disorders_at_high_T(tmp_path):
    """2D XY ferromagnet (BKT scale T_BKT ≈ 0.9): on a small lattice the
    magnetisation is high deep below it and O(N^{-1/2}) far above it."""

    def tail_magn(T):
        vals = []
        for sd in (1, 2, 3):
            sg = make_lattice(tmp_path, pflip=0.0)
            m = XYMetropolis(sg, T=T, steps=200, seed=sd, savedisk=False)
            m.run()
            vals.append(float(np.mean(m.magn[-40:])))
        return float(np.mean(vals))

    assert tail_magn(0.2) > 0.5
    assert tail_magn(5.0) < 0.35


def test_sync_mode_checkerboard_and_equilibrium_agreement(tmp_path):
    """The shared greedy coloring gives the 2-class checkerboard on the
    bipartite square lattice, and sync agrees with async in equilibrium."""
    sg = make_lattice(tmp_path, pflip=0.0)
    m_sync = XYMetropolis(
        sg, T=2.0, upd_mode="sync", steps=300, seed=SEED, savedisk=False
    )
    m_sync.run()
    assert len(m_sync._engine.color_classes) == 2
    m_async = XYMetropolis(sg, T=2.0, steps=300, seed=SEED, savedisk=False)
    m_async.run()
    e_sync = float(np.mean(m_sync.ene[150:]))
    e_async = float(np.mean(m_async.ene[150:]))
    assert abs(e_sync - e_async) < 0.1


def test_exchange_move_conserves_angle_multiset(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = XYMetropolis(
        sg, T=1.0, move="exchange", steps=10, seed=SEED, savedisk=False
    )
    m.init_dynamics()
    angles_before = np.sort(m.s.copy())
    m.run()
    assert np.allclose(np.sort(m.s), angles_before)
    assert np.isclose(m._e_running, m.total_energy())


def test_wolff_orders_ferromagnet(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m = XYMetropolis(
        sg, T=0.3, move="wolff", steps=20, seed=SEED, savedisk=False
    )
    m.run()
    assert float(np.mean(m.magn[-5:])) > 0.6
    assert m.cluster_sizes is not None and len(m.cluster_sizes) == 20
    assert np.isclose(m._e_running, m.total_energy())


def test_sw_runs_ferromagnet(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m = XYMetropolis(sg, T=1.0, move="sw", steps=15, seed=SEED, savedisk=False)
    m.run()
    assert m.cluster_counts is not None and len(m.cluster_counts) == 15
    assert np.isclose(m._e_running, m.total_energy())


def test_cluster_moves_valid_on_signed_couplings(tmp_path):
    """Unlike Potts, the XY cluster moves accept signed graphs: the
    embedded-Ising reflection (ref M6) is an involution, so the bond
    satisfaction is flip-covariant exactly as in the Ising case. Pin that
    both engines RUN on a signed graph and keep the energy bookkeeping."""
    for move in ("wolff", "sw"):
        sg = make_lattice(tmp_path, pflip=0.15)
        m = XYMetropolis(
            sg, T=0.8, move=move, steps=8, seed=SEED, savedisk=False
        )
        m.run()
        assert np.isclose(m._e_running, m.total_energy())


# =========================================================================
# Capability boundaries (hard errors, invariant #3)
# =========================================================================
def test_continuous_state_axis_exclusions(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(NotImplementedError, match="von-Mises"):
        XYMetropolis(sg, rule="heatbath", savedisk=False)
    with pytest.raises(NotImplementedError, match="deterministic local"):
        XYMetropolis(sg, upd_mode="gillespie", savedisk=False)
    with pytest.raises(ValueError, match="measure zero"):
        XYMetropolis(sg, tie_flip_p=0.5, savedisk=False)
    with pytest.raises(ValueError, match="measure zero"):
        XYMetropolis(sg, tie_flip_p="frozen", savedisk=False)


def test_axis_capability_errors(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="delta must be > 0"):
        XYMetropolis(sg, delta=0.0, savedisk=False)
    with pytest.raises(ValueError, match="rule"):
        XYMetropolis(sg, rule="bogus", savedisk=False)
    with pytest.raises(ValueError, match="move"):
        XYMetropolis(sg, move="bogus", savedisk=False)
    # No spectral move for XY: the eigenvector-subspace walk is an
    # Ising (±1) construction.
    with pytest.raises(ValueError, match="move"):
        XYMetropolis(sg, move="spectral", savedisk=False)
    with pytest.raises(ValueError, match="order does not apply"):
        XYMetropolis(sg, move="wolff", order="typewriter", savedisk=False)
    with pytest.raises(ValueError, match="rule does not apply"):
        XYMetropolis(sg, move="sw", rule="glauber", savedisk=False)
    with pytest.raises(ValueError, match="delta .* does not apply"):
        XYMetropolis(sg, move="wolff", delta=0.7, savedisk=False)
    with pytest.raises(ValueError, match="asynchronous kinetics"):
        XYMetropolis(sg, move="exchange", upd_mode="sync", savedisk=False)
    with pytest.raises(ValueError, match="observable"):
        XYMetropolis(sg, observables=("bogus",), savedisk=False)
    with pytest.raises(ValueError, match="coupling_norm"):
        XYMetropolis(sg, coupling_norm="bogus", savedisk=False)
    with pytest.raises(NotImplementedError, match="field"):
        XYMetropolis(sg, field=np.ones(sg.N), savedisk=False)


def test_pb_capability_guards(tmp_path):
    """The pb solver's supports() gate rejects, at setup, everything the
    compiled kernel cannot faithfully represent (no binary needed)."""
    sg = make_lattice(tmp_path)

    def pb(**kw):
        m = XYMetropolis(
            sg, T=1.0, steps=2, seed=SEED, savedisk=False, runlang="pb", **kw
        )
        m.run()

    with pytest.raises(NotImplementedError, match="move='single'"):
        pb(move="exchange")
    with pytest.raises(NotImplementedError, match="asynchronous only"):
        pb(upd_mode="sync")
    with pytest.raises(NotImplementedError, match="metropolis"):
        pb(rule="glauber")
    with pytest.raises(NotImplementedError, match="random"):
        pb(order="typewriter")
    with pytest.raises(NotImplementedError, match="snapshots"):
        pb(observables=("energy", "magn", "snapshots"))


# =========================================================================
# Native pybind backend (kernel semantics; skipped without the binary)
# =========================================================================
@needs_native
def test_pb_bookkeeping_and_reproducibility(tmp_path):
    sg = make_lattice(tmp_path)
    runs = []
    for _ in range(2):
        m = XYMetropolis(
            sg, T=1.0, steps=6, seed=SEED, savedisk=False, runlang="pb"
        )
        m.run()
        # exact record-then-sweep reconstruction + running-energy invariant
        assert np.isclose(m.ene[0], m.energy_intensive(m.sini))
        assert np.isclose(m._e_running, m.total_energy())
        runs.append((m.s.copy(), np.asarray(m.ene).copy()))
    assert np.array_equal(runs[0][0], runs[1][0])
    assert np.allclose(runs[0][1], runs[1][1])
    other = XYMetropolis(
        sg, T=1.0, steps=6, seed=SEED + 1, savedisk=False, runlang="pb"
    )
    other.run()
    assert not np.array_equal(runs[0][0], other.s)


@needs_native
def test_pb_zero_T_is_monotone(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = XYMetropolis(
        sg, T=0.0, steps=20, seed=SEED, savedisk=False, runlang="pb"
    )
    m.run()
    assert np.all(np.diff(np.asarray(m.ene)) <= 1e-12)


@needs_native
def test_pb_matches_py_order_disorder(tmp_path):
    """Different RNG streams — compare the qualitative phase (D-B6:
    statistical, never byte-equality)."""

    def tail_magn(runlang, T):
        vals = []
        for sd in (1, 2, 3):
            sg = make_lattice(tmp_path, pflip=0.0)
            m = XYMetropolis(
                sg, T=T, steps=200, seed=sd, savedisk=False, runlang=runlang
            )
            m.run()
            vals.append(float(np.mean(m.magn[-40:])))
        return float(np.mean(vals))

    assert tail_magn("pb", 0.2) > 0.5
    assert tail_magn("pb", 5.0) < 0.35
