"""Heisenberg vertical slice: HeisenbergBase + HeisenbergMetropolis.

The fourth consumer of the statsys-unification architecture (D-B8/D-B9,
after Ising/Potts/XY) and the first VECTOR-state substrate (a unit 3-vector
per site): the SAME ``statsys._thermal`` engines drive an (N, 3) float64
state — the only engine change was making the Kawasaki identity-swap check
state-shape-agnostic. Pins:

- the §3b Hamiltonian bookkeeping (ΔE↔E from ONE symmetric J; incremental
  running energy == full recompute after every run);
- the exact Ising limit (spins restricted to ±z give ``n_i · n_j = s_i s_j``,
  so Heisenberg energies and z-flip ΔE's equal the Ising ones);
- physics sanity per move/update mode (quench monotonicity, order/disorder,
  spin-multiset conservation under Kawasaki, cluster ordering — including on
  SIGNED graphs via the embedded-Ising reflection, ref M6);
- the capability boundaries (heatbath/gillespie not wired for a continuous
  state, tie_flip_p meaningless, pb kernel guards) — hard errors,
  invariant #3;
- the legacy HeisenbergModel stays untouched and agrees on the Hamiltonian.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.HeisenbergModel import (
    HeisenbergMetropolis,
    HeisenbergModel,
)
from lrgsglib.statsys.IsingDynamics import IsingMetropolis

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
        from lrgsglib.statsys.HeisenbergModel.ccore import (  # noqa: F401
            _heisenberg_native,
        )

        return True
    except Exception:
        return False


needs_native = pytest.mark.skipif(
    not _native_available(), reason="_heisenberg_native not built"
)


# =========================================================================
# Bookkeeping: the §3b contract on the vector substrate
# =========================================================================
def test_defaults_construct_and_run(tmp_path):
    sg = make_lattice(tmp_path)
    model = HeisenbergMetropolis(sg, T=1.0, steps=5, seed=SEED, savedisk=False)
    model.run()
    assert len(model.ene) == 5 and len(model.magn) == 5
    assert model.s.shape == (N, 3) and model.s.dtype == np.float64
    assert np.allclose(np.linalg.norm(model.s, axis=1), 1.0)
    # record-then-sweep cadence + incremental energy invariant
    assert np.isclose(model.ene[0], model.energy_intensive(model.sini))
    assert np.isclose(model._e_running, model.total_energy())


def test_hamiltonian_matches_legacy(tmp_path):
    """The new CSR total_energy equals the legacy dense-loop energy()."""
    sg = make_lattice(tmp_path, pflip=0.2)
    new = HeisenbergMetropolis(sg, T=1.0, steps=1, seed=SEED, savedisk=False)
    new.init_dynamics()
    legacy = HeisenbergModel(sg, T=1.0, steps=1, seed=SEED)
    assert np.isclose(new.total_energy(), legacy.energy(new.s))


def test_axis_spins_are_ising(tmp_path):
    """Spins restricted to ±z give n_i · n_j = s_i s_j, so the Heisenberg
    Hamiltonian and the z-flip ΔE equal the Ising ones on the same signed
    graph, exactly."""
    sg = make_lattice(tmp_path, pflip=0.15)
    rng = np.random.RandomState(3)
    spins = rng.choice([-1, 1], size=sg.N).astype(np.int8)
    vecs = np.zeros((sg.N, 3))
    vecs[:, 2] = spins
    hb = HeisenbergMetropolis(
        sg, T=1.0, steps=1, seed=3, ic="custom", savedisk=False
    )
    hb.init_dynamics(custom=vecs)
    ising = IsingMetropolis(sg, T=1.0, steps=1, seed=3, savedisk=False)
    ising.init_dynamics()
    ising.s = spins
    assert np.isclose(hb.total_energy(), ising.total_energy())
    for nd in range(sg.N):
        dE_hb = hb.delta_E(nd, -vecs[nd])
        dE_is = ising.delta_E(nd, -int(spins[nd]))
        assert np.isclose(dE_hb, dE_is)


def test_seeded_reproducibility(tmp_path):
    finals = []
    for _ in range(2):
        sg = make_lattice(tmp_path)
        m = HeisenbergMetropolis(sg, T=1.5, steps=6, seed=SEED, savedisk=False)
        m.run()
        finals.append(m.s.copy())
    assert np.array_equal(finals[0], finals[1])


# =========================================================================
# Physics: quench, order/disorder, update modes
# =========================================================================
def test_zero_T_quench_is_monotone(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = HeisenbergMetropolis(sg, T=0.0, steps=20, seed=SEED, savedisk=False)
    m.run()
    assert np.all(np.diff(np.asarray(m.ene)) <= 1e-12)


def test_orders_at_low_T_disorders_at_high_T(tmp_path):
    """On a small 2D lattice deep below the coupling scale the correlation
    length exceeds the box (finite-size order despite Mermin-Wagner); far
    above it the magnetisation is O(N^{-1/2})."""

    def tail_magn(T):
        vals = []
        for sd in (1, 2, 3):
            sg = make_lattice(tmp_path, pflip=0.0)
            m = HeisenbergMetropolis(
                sg, T=T, steps=200, seed=sd, savedisk=False
            )
            m.run()
            vals.append(float(np.mean(m.magn[-40:])))
        return float(np.mean(vals))

    assert tail_magn(0.1) > 0.5
    assert tail_magn(5.0) < 0.35


def test_sync_mode_checkerboard_and_equilibrium_agreement(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m_sync = HeisenbergMetropolis(
        sg, T=2.0, upd_mode="sync", steps=300, seed=SEED, savedisk=False
    )
    m_sync.run()
    assert len(m_sync._engine.color_classes) == 2
    m_async = HeisenbergMetropolis(
        sg, T=2.0, steps=300, seed=SEED, savedisk=False
    )
    m_async.run()
    e_sync = float(np.mean(m_sync.ene[150:]))
    e_async = float(np.mean(m_async.ene[150:]))
    assert abs(e_sync - e_async) < 0.1


def test_exchange_move_conserves_spin_multiset(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = HeisenbergMetropolis(
        sg, T=1.0, move="exchange", steps=10, seed=SEED, savedisk=False
    )
    m.init_dynamics()
    # The swap permutes rows: compare the row multisets via lexsort.
    before = m.s.copy()
    m.run()
    order_b = np.lexsort(before.T)
    order_a = np.lexsort(m.s.T)
    assert np.allclose(before[order_b], m.s[order_a])
    assert np.isclose(m._e_running, m.total_energy())


def test_wolff_orders_ferromagnet(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m = HeisenbergMetropolis(
        sg, T=0.2, move="wolff", steps=25, seed=SEED, savedisk=False
    )
    m.run()
    assert float(np.mean(m.magn[-5:])) > 0.6
    assert m.cluster_sizes is not None and len(m.cluster_sizes) == 25
    assert np.isclose(m._e_running, m.total_energy())


def test_sw_runs_ferromagnet(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m = HeisenbergMetropolis(
        sg, T=1.0, move="sw", steps=15, seed=SEED, savedisk=False
    )
    m.run()
    assert m.cluster_counts is not None and len(m.cluster_counts) == 15
    assert np.isclose(m._e_running, m.total_energy())


def test_cluster_moves_valid_on_signed_couplings(tmp_path):
    """The embedded-Ising reflection (ref M6) is an involution, so the bond
    satisfaction is flip-covariant and cluster moves stay valid on signed
    graphs. Pin that both engines RUN and keep the energy bookkeeping."""
    for move in ("wolff", "sw"):
        sg = make_lattice(tmp_path, pflip=0.15)
        m = HeisenbergMetropolis(
            sg, T=0.8, move=move, steps=8, seed=SEED, savedisk=False
        )
        m.run()
        assert np.isclose(m._e_running, m.total_energy())


# =========================================================================
# Capability boundaries (hard errors, invariant #3)
# =========================================================================
def test_continuous_state_axis_exclusions(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(NotImplementedError, match="conditional resample"):
        HeisenbergMetropolis(sg, rule="heatbath", savedisk=False)
    with pytest.raises(NotImplementedError, match="deterministic local"):
        HeisenbergMetropolis(sg, upd_mode="gillespie", savedisk=False)
    with pytest.raises(ValueError, match="measure zero"):
        HeisenbergMetropolis(sg, tie_flip_p=0.5, savedisk=False)


def test_axis_capability_errors(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="delta must be > 0"):
        HeisenbergMetropolis(sg, delta=0.0, savedisk=False)
    with pytest.raises(ValueError, match="rule"):
        HeisenbergMetropolis(sg, rule="bogus", savedisk=False)
    with pytest.raises(ValueError, match="move"):
        HeisenbergMetropolis(sg, move="bogus", savedisk=False)
    with pytest.raises(ValueError, match="move"):
        HeisenbergMetropolis(sg, move="spectral", savedisk=False)
    with pytest.raises(ValueError, match="order does not apply"):
        HeisenbergMetropolis(
            sg, move="wolff", order="typewriter", savedisk=False
        )
    with pytest.raises(ValueError, match="rule does not apply"):
        HeisenbergMetropolis(sg, move="sw", rule="glauber", savedisk=False)
    with pytest.raises(ValueError, match="delta .* does not apply"):
        HeisenbergMetropolis(sg, move="wolff", delta=0.7, savedisk=False)
    with pytest.raises(ValueError, match="asynchronous kinetics"):
        HeisenbergMetropolis(
            sg, move="exchange", upd_mode="sync", savedisk=False
        )
    with pytest.raises(ValueError, match="observable"):
        HeisenbergMetropolis(sg, observables=("bogus",), savedisk=False)
    with pytest.raises(ValueError, match="coupling_norm"):
        HeisenbergMetropolis(sg, coupling_norm="bogus", savedisk=False)
    with pytest.raises(NotImplementedError, match="field"):
        HeisenbergMetropolis(sg, field=np.ones(sg.N), savedisk=False)


def test_pb_capability_guards(tmp_path):
    """The pb solver's supports() gate rejects, at setup, everything the
    compiled kernel cannot faithfully represent (no binary needed)."""
    sg = make_lattice(tmp_path)

    def pb(**kw):
        m = HeisenbergMetropolis(
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
        m = HeisenbergMetropolis(
            sg, T=1.0, steps=6, seed=SEED, savedisk=False, runlang="pb"
        )
        m.run()
        # exact record-then-sweep reconstruction + running-energy invariant
        assert np.isclose(m.ene[0], m.energy_intensive(m.sini))
        assert np.isclose(m._e_running, m.total_energy())
        runs.append((m.s.copy(), np.asarray(m.ene).copy()))
    assert np.array_equal(runs[0][0], runs[1][0])
    assert np.allclose(runs[0][1], runs[1][1])
    other = HeisenbergMetropolis(
        sg, T=1.0, steps=6, seed=SEED + 1, savedisk=False, runlang="pb"
    )
    other.run()
    assert not np.array_equal(runs[0][0], other.s)


@needs_native
def test_pb_zero_T_is_monotone(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = HeisenbergMetropolis(
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
            m = HeisenbergMetropolis(
                sg, T=T, steps=200, seed=sd, savedisk=False, runlang=runlang
            )
            m.run()
            vals.append(float(np.mean(m.magn[-40:])))
        return float(np.mean(vals))

    assert tail_magn("pb", 0.1) > 0.5
    assert tail_magn("pb", 5.0) < 0.35
