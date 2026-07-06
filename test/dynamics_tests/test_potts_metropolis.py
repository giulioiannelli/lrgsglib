"""Potts vertical slice: PottsBase + PottsMetropolis on the shared spine.

The second consumer of the statsys-unification architecture (D-B8/D-B9,
Potts after Ising): the SAME ``statsys._thermal`` engines drive a different
substrate — no engine code changed. Pins:

- the §3b Hamiltonian bookkeeping (ΔE↔E from ONE symmetric J; incremental
  running energy == full recompute after every run);
- the exact q = 2 ↔ Ising correspondence (``δ = (1 + s_i s_j)/2`` → Potts
  energies are Ising energies at halved couplings, up to a constant);
- physics sanity per move/update mode (quench monotonicity, order/disorder,
  density conservation under Kawasaki, cluster ordering on ferro graphs);
- the capability boundaries (heatbath/gillespie q = 2 only, cluster moves
  non-negative couplings only, pb kernel guards) — hard errors, invariant #3;
- the legacy PottsModel stays untouched and agrees on the Hamiltonian.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import IsingMetropolis
from lrgsglib.statsys.PottsModel import PottsMetropolis, PottsModel

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
        from lrgsglib.statsys.PottsModel.ccore import (  # noqa: F401
            _potts_native,
        )

        return True
    except Exception:
        return False


needs_native = pytest.mark.skipif(
    not _native_available(), reason="_potts_native not built"
)


# =========================================================================
# Bookkeeping: the §3b contract on the new substrate
# =========================================================================
def test_defaults_construct_and_run(tmp_path):
    sg = make_lattice(tmp_path)
    model = PottsMetropolis(sg, q=3, T=1.0, steps=5, seed=SEED, savedisk=False)
    model.run()
    assert len(model.ene) == 5 and len(model.magn) == 5
    assert model.s.dtype == np.int32
    assert ((model.s >= 0) & (model.s < 3)).all()
    # record-then-sweep cadence + incremental energy invariant
    assert np.isclose(model.ene[0], model.energy_intensive(model.sini))
    assert np.isclose(model._e_running, model.total_energy())


def test_hamiltonian_matches_legacy(tmp_path):
    """The new CSR total_energy equals the legacy dense-loop energy()."""
    sg = make_lattice(tmp_path, pflip=0.2)
    new = PottsMetropolis(sg, q=3, T=1.0, steps=1, seed=SEED, savedisk=False)
    new.init_dynamics()
    legacy = PottsModel(sg, q=3, T=1.0, steps=1, seed=SEED)
    assert np.isclose(new.total_energy(), legacy.energy(new.s))


def test_q2_is_ising_at_halved_couplings(tmp_path):
    """δ(σi,σj) = (1 + s_i s_j)/2 under σ → s = 2σ − 1 (ref M10):
    E_potts = E_ising/2 − ΣJ/2 and ΔE_potts = ΔE_ising/2, exactly."""
    sg = make_lattice(tmp_path, pflip=0.15)
    potts = PottsMetropolis(sg, q=2, T=1.0, steps=1, seed=3, savedisk=False)
    potts.init_dynamics()
    ising = IsingMetropolis(sg, T=2.0, steps=1, seed=3, savedisk=False)
    ising.init_dynamics()
    ising.s = (2 * potts.s - 1).astype(np.int8)
    sum_J = potts._edge_J.sum()
    assert np.isclose(
        potts.total_energy(), ising.total_energy() / 2.0 - sum_J / 2.0
    )
    for nd in range(sg.N):
        dE_p = potts.delta_E(nd, 1 - int(potts.s[nd]))
        dE_i = ising.delta_E(nd, -int(ising.s[nd]))
        assert np.isclose(dE_p, dE_i / 2.0)


def test_seeded_reproducibility(tmp_path):
    finals = []
    for _ in range(2):
        sg = make_lattice(tmp_path)
        m = PottsMetropolis(sg, q=3, T=1.5, steps=6, seed=SEED, savedisk=False)
        m.run()
        finals.append(m.s.copy())
    assert np.array_equal(finals[0], finals[1])


# =========================================================================
# Physics: quench, order/disorder, update modes
# =========================================================================
def test_zero_T_quench_is_monotone_with_frozen_ties(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = PottsMetropolis(
        sg,
        q=3,
        T=0.0,
        tie_flip_p="frozen",
        steps=20,
        seed=SEED,
        savedisk=False,
    )
    m.run()
    assert np.all(np.diff(np.asarray(m.ene)) <= 1e-12)


def test_orders_at_low_T_disorders_at_high_T(tmp_path):
    """q = 3 ferromagnet (2D Tc = 1/ln(1+√3) ≈ 0.995, ref M10): the order
    parameter is high deep below Tc and small deep above it."""

    def tail_order(T):
        vals = []
        for sd in (1, 2, 3):
            sg = make_lattice(tmp_path, pflip=0.0)
            m = PottsMetropolis(
                sg, q=3, T=T, steps=150, seed=sd, savedisk=False
            )
            m.run()
            vals.append(float(np.mean(m.magn[-30:])))
        return float(np.mean(vals))

    assert tail_order(0.2) > 0.5
    assert tail_order(5.0) < 0.3


def test_sync_mode_checkerboard_and_equilibrium_agreement(tmp_path):
    """The shared greedy coloring gives the 2-class checkerboard on the
    bipartite square lattice, and sync agrees with async in equilibrium."""
    sg = make_lattice(tmp_path, pflip=0.0)
    m_sync = PottsMetropolis(
        sg, q=3, T=2.0, upd_mode="sync", steps=300, seed=SEED, savedisk=False
    )
    m_sync.run()
    assert len(m_sync._engine.color_classes) == 2
    m_async = PottsMetropolis(
        sg, q=3, T=2.0, steps=300, seed=SEED, savedisk=False
    )
    m_async.run()
    e_sync = float(np.mean(m_sync.ene[150:]))
    e_async = float(np.mean(m_async.ene[150:]))
    assert abs(e_sync - e_async) < 0.1


def test_gillespie_q2_freezes_at_zero_T(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = PottsMetropolis(
        sg,
        q=2,
        T=0.0,
        tie_flip_p="frozen",
        upd_mode="gillespie",
        steps=30,
        seed=SEED,
        savedisk=False,
    )
    m.run()
    assert m.is_absorbing()
    ene = np.asarray(m.ene)
    assert np.all(np.diff(ene) <= 1e-12)


def test_exchange_move_conserves_state_densities(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = PottsMetropolis(
        sg, q=3, T=1.0, move="exchange", steps=10, seed=SEED, savedisk=False
    )
    m.init_dynamics()
    counts_before = np.bincount(m.s, minlength=3).copy()
    m.run()
    assert np.array_equal(np.bincount(m.s, minlength=3), counts_before)
    assert np.isclose(m._e_running, m.total_energy())


def test_wolff_orders_ferromagnet(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m = PottsMetropolis(
        sg, q=3, T=0.5, move="wolff", steps=20, seed=SEED, savedisk=False
    )
    m.run()
    assert float(np.mean(m.magn[-5:])) > 0.6
    assert m.cluster_sizes is not None and len(m.cluster_sizes) == 20
    assert np.isclose(m._e_running, m.total_energy())


def test_sw_runs_ferromagnet(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m = PottsMetropolis(
        sg, q=3, T=1.5, move="sw", steps=15, seed=SEED, savedisk=False
    )
    m.run()
    assert m.cluster_counts is not None and len(m.cluster_counts) == 15
    assert np.isclose(m._e_running, m.total_energy())


# =========================================================================
# Capability boundaries (hard errors, invariant #3)
# =========================================================================
def test_cluster_moves_reject_signed_couplings(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    for move in ("wolff", "sw"):
        m = PottsMetropolis(
            sg, q=3, T=1.0, move=move, steps=2, seed=SEED, savedisk=False
        )
        with pytest.raises(NotImplementedError, match="non-negative"):
            m.run()


def test_two_state_only_axes(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(NotImplementedError, match="two-state"):
        PottsMetropolis(sg, q=3, rule="heatbath", savedisk=False)
    with pytest.raises(NotImplementedError, match="two-state"):
        PottsMetropolis(sg, q=3, upd_mode="gillespie", savedisk=False)
    # ... and both are available at q = 2.
    PottsMetropolis(sg, q=2, rule="heatbath", savedisk=False)
    PottsMetropolis(sg, q=2, upd_mode="gillespie", savedisk=False)


def test_axis_capability_errors(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="q must be >= 2"):
        PottsMetropolis(sg, q=1, savedisk=False)
    with pytest.raises(ValueError, match="rule"):
        PottsMetropolis(sg, rule="bogus", savedisk=False)
    with pytest.raises(ValueError, match="move"):
        PottsMetropolis(sg, move="bogus", savedisk=False)
    # No spectral move for Potts: the eigenvector-subspace walk is an
    # Ising (±1) construction.
    with pytest.raises(ValueError, match="move"):
        PottsMetropolis(sg, move="spectral", savedisk=False)
    with pytest.raises(ValueError, match="order applies"):
        PottsMetropolis(
            sg, q=2, upd_mode="gillespie", order="typewriter", savedisk=False
        )
    with pytest.raises(ValueError, match="order does not apply"):
        PottsMetropolis(sg, move="wolff", order="typewriter", savedisk=False)
    with pytest.raises(ValueError, match="rule does not apply"):
        PottsMetropolis(sg, move="sw", rule="glauber", savedisk=False)
    with pytest.raises(ValueError, match="tie_flip_p"):
        PottsMetropolis(sg, move="wolff", tie_flip_p=0.5, savedisk=False)
    with pytest.raises(ValueError, match="asynchronous kinetics"):
        PottsMetropolis(sg, move="exchange", upd_mode="sync", savedisk=False)
    with pytest.raises(ValueError, match="tie_flip_p is metropolis-only"):
        PottsMetropolis(sg, rule="glauber", tie_flip_p=0.5, savedisk=False)
    with pytest.raises(ValueError, match="observable"):
        PottsMetropolis(sg, observables=("bogus",), savedisk=False)
    with pytest.raises(ValueError, match="coupling_norm"):
        PottsMetropolis(sg, coupling_norm="bogus", savedisk=False)
    with pytest.raises(NotImplementedError, match="field"):
        PottsMetropolis(sg, field=np.ones(sg.N), savedisk=False)


def test_pb_capability_guards(tmp_path):
    """The pb solver's supports() gate rejects, at setup, everything the
    compiled kernel cannot faithfully represent (no binary needed)."""
    sg = make_lattice(tmp_path)

    def pb(**kw):
        m = PottsMetropolis(
            sg,
            q=3,
            T=1.0,
            steps=2,
            seed=SEED,
            savedisk=False,
            runlang="pb",
            **kw,
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
    with pytest.raises(NotImplementedError, match="tie"):
        pb(tie_flip_p=0.5)
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
        m = PottsMetropolis(
            sg, q=3, T=1.0, steps=6, seed=SEED, savedisk=False, runlang="pb"
        )
        m.run()
        # exact record-then-sweep reconstruction + running-energy invariant
        assert np.isclose(m.ene[0], m.energy_intensive(m.sini))
        assert np.isclose(m._e_running, m.total_energy())
        runs.append((m.s.copy(), np.asarray(m.ene).copy()))
    assert np.array_equal(runs[0][0], runs[1][0])
    assert np.allclose(runs[0][1], runs[1][1])
    other = PottsMetropolis(
        sg, q=3, T=1.0, steps=6, seed=SEED + 1, savedisk=False, runlang="pb"
    )
    other.run()
    assert not np.array_equal(runs[0][0], other.s)


@needs_native
def test_pb_zero_T_is_monotone(tmp_path):
    """The kernel accepts every ΔE <= 0 proposal, so at T = 0 the energy
    trace is nonincreasing (ties flip freely but cost nothing)."""
    sg = make_lattice(tmp_path, pflip=0.1)
    m = PottsMetropolis(
        sg, q=3, T=0.0, steps=20, seed=SEED, savedisk=False, runlang="pb"
    )
    m.run()
    assert np.all(np.diff(np.asarray(m.ene)) <= 1e-12)


@needs_native
def test_pb_matches_py_order_disorder(tmp_path):
    """Different RNG streams — compare the qualitative phase (D-B6:
    statistical, never byte-equality)."""

    def tail_order(runlang, T):
        vals = []
        for sd in (1, 2, 3):
            sg = make_lattice(tmp_path, pflip=0.0)
            m = PottsMetropolis(
                sg,
                q=3,
                T=T,
                steps=150,
                seed=sd,
                savedisk=False,
                runlang=runlang,
            )
            m.run()
            vals.append(float(np.mean(m.magn[-30:])))
        return float(np.mean(vals))

    assert tail_order("pb", 0.2) > 0.5
    assert tail_order("pb", 5.0) < 0.3
