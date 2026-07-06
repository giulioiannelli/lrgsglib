"""MultiSpecies vertical slice: MultiSpeciesBase + MultiSpeciesMetropolis.

The last member of the spin quad (statsys unification D-B8/D-B9): k coupled
Potts components per site on the SAME ``statsys._thermal`` engines — no
engine code changed. Pins:

- the §3b Hamiltonian bookkeeping (ΔE↔E from ONE symmetric J and ONE
  symmetric interaction matrix; incremental running energy == full recompute
  after every run);
- the exact Potts limit (k = 1, M = 1 → q-state Potts, energies and ΔE's
  equal PottsMetropolis, ref M10);
- physics sanity per move/update mode (quench monotonicity, order/disorder,
  per-species density conservation under the full-site Kawasaki swap);
- the capability boundaries (asymmetric interaction matrices, cluster
  moves, heatbath/gillespie, pb kernel restrictions — identity M, uniform
  q, energy-only records) — hard errors, invariant #3;
- the legacy MultiSpeciesModel stays untouched and agrees on the
  Hamiltonian (identity AND symmetric non-identity M).
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.MultiSpeciesModel import (
    MultiSpeciesMetropolis,
    MultiSpeciesModel,
)
from lrgsglib.statsys.PottsModel import PottsMetropolis

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
        from lrgsglib.statsys.MultiSpeciesModel.ccore import (  # noqa: F401
            _multispec_native,
        )

        return True
    except Exception:
        return False


needs_native = pytest.mark.skipif(
    not _native_available(), reason="_multispec_native not built"
)


# =========================================================================
# Bookkeeping: the §3b contract on the k-component substrate
# =========================================================================
def test_defaults_construct_and_run(tmp_path):
    sg = make_lattice(tmp_path)
    model = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=1.0,
        steps=5,
        seed=SEED,
        savedisk=False,
    )
    model.run()
    assert len(model.ene) == 5 and len(model.magn) == 5
    assert model.s.shape == (N, 2) and model.s.dtype == np.int32
    assert ((model.s >= 0) & (model.s < 3)).all()
    # record-then-sweep cadence + incremental energy invariant
    assert np.isclose(model.ene[0], model.energy_intensive(model.sini))
    assert np.isclose(model._e_running, model.total_energy())


def test_hamiltonian_matches_legacy_identity_M(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.2)
    new = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=1.0,
        steps=1,
        seed=SEED,
        savedisk=False,
    )
    new.init_dynamics()
    legacy = MultiSpeciesModel(sg, species=2, q_per_species=3, T=1.0, steps=1)
    assert np.isclose(new.total_energy(), legacy.energy(new.s))


def test_hamiltonian_matches_legacy_symmetric_M(tmp_path):
    """A symmetric non-identity interaction matrix: the new vectorized
    energy equals the legacy quadruple loop."""
    sg = make_lattice(tmp_path, pflip=0.2)
    M = np.array([[1.0, 0.4], [0.4, 0.7]])
    new = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=1.0,
        interaction_matrix=M,
        steps=1,
        seed=SEED,
        savedisk=False,
    )
    new.init_dynamics()
    legacy = MultiSpeciesModel(
        sg, species=2, q_per_species=3, T=1.0, interaction_matrix=M, steps=1
    )
    assert np.isclose(new.total_energy(), legacy.energy(new.s))


def test_single_species_is_potts(tmp_path):
    """k = 1 with M = 1 is exactly the q-state Potts model (ref M10):
    same energy and same per-site ΔE as PottsMetropolis."""
    sg = make_lattice(tmp_path, pflip=0.15)
    ms = MultiSpeciesMetropolis(
        sg,
        species=1,
        q_per_species=3,
        T=1.0,
        steps=1,
        seed=3,
        savedisk=False,
    )
    ms.init_dynamics()
    potts = PottsMetropolis(sg, q=3, T=1.0, steps=1, seed=3, savedisk=False)
    potts.init_dynamics()
    potts.s = ms.s[:, 0].copy()
    assert np.isclose(ms.total_energy(), potts.total_energy())
    for nd in range(sg.N):
        cur = int(ms.s[nd, 0])
        for label in range(3):
            if label == cur:
                continue
            assert np.isclose(
                ms.delta_E(nd, (0, label)), potts.delta_E(nd, label)
            )


def test_seeded_reproducibility(tmp_path):
    finals = []
    for _ in range(2):
        sg = make_lattice(tmp_path)
        m = MultiSpeciesMetropolis(
            sg,
            species=2,
            q_per_species=3,
            T=1.5,
            steps=6,
            seed=SEED,
            savedisk=False,
        )
        m.run()
        finals.append(m.s.copy())
    assert np.array_equal(finals[0], finals[1])


# =========================================================================
# Physics: quench, order/disorder, update modes
# =========================================================================
def test_zero_T_quench_is_monotone_with_frozen_ties(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=0.0,
        tie_flip_p="frozen",
        steps=20,
        seed=SEED,
        savedisk=False,
    )
    m.run()
    assert np.all(np.diff(np.asarray(m.ene)) <= 1e-12)


def test_orders_at_low_T_disorders_at_high_T(tmp_path):
    """Identity M: each species is an independent q = 3 Potts layer
    (2D Tc ≈ 0.995, ref M10) — the species-averaged order parameter is
    high deep below Tc and small deep above it."""

    def tail_order(T):
        vals = []
        for sd in (1, 2, 3):
            sg = make_lattice(tmp_path, pflip=0.0)
            m = MultiSpeciesMetropolis(
                sg,
                species=2,
                q_per_species=3,
                T=T,
                steps=150,
                seed=sd,
                savedisk=False,
            )
            m.run()
            vals.append(float(np.mean(m.magn[-30:])))
        return float(np.mean(vals))

    assert tail_order(0.2) > 0.5
    assert tail_order(5.0) < 0.3


def test_sync_mode_checkerboard_and_equilibrium_agreement(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.0)
    m_sync = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=2.0,
        upd_mode="sync",
        steps=300,
        seed=SEED,
        savedisk=False,
    )
    m_sync.run()
    assert len(m_sync._engine.color_classes) == 2
    m_async = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=2.0,
        steps=300,
        seed=SEED,
        savedisk=False,
    )
    m_async.run()
    e_sync = float(np.mean(m_sync.ene[150:]))
    e_async = float(np.mean(m_async.ene[150:]))
    assert abs(e_sync - e_async) < 0.15


def test_exchange_move_conserves_per_species_densities(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    m = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=(3, 4),
        T=1.0,
        move="exchange",
        steps=10,
        seed=SEED,
        savedisk=False,
    )
    m.init_dynamics()
    before = [np.bincount(m.s[:, k], minlength=5).copy() for k in range(2)]
    m.run()
    for k in range(2):
        assert np.array_equal(np.bincount(m.s[:, k], minlength=5), before[k])
    assert np.isclose(m._e_running, m.total_energy())


def test_per_species_q_heterogeneous(tmp_path):
    """Different q per species: labels stay in each species' own range and
    the bookkeeping holds."""
    sg = make_lattice(tmp_path, pflip=0.1)
    m = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=(2, 5),
        T=1.0,
        steps=10,
        seed=SEED,
        savedisk=False,
    )
    m.run()
    assert ((m.s[:, 0] >= 0) & (m.s[:, 0] < 2)).all()
    assert ((m.s[:, 1] >= 0) & (m.s[:, 1] < 5)).all()
    assert np.isclose(m._e_running, m.total_energy())


def test_coupled_species_interaction_runs(tmp_path):
    """A symmetric off-diagonal M couples the species; the incremental
    energy stays consistent with the full Hamiltonian."""
    sg = make_lattice(tmp_path, pflip=0.1)
    M = np.array([[1.0, 0.5], [0.5, 1.0]])
    m = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=1.0,
        interaction_matrix=M,
        steps=15,
        seed=SEED,
        savedisk=False,
    )
    m.run()
    assert np.isclose(m._e_running, m.total_energy())


# =========================================================================
# Capability boundaries (hard errors, invariant #3)
# =========================================================================
def test_asymmetric_interaction_matrix_rejected(tmp_path):
    sg = make_lattice(tmp_path)
    M = np.array([[1.0, 0.8], [0.2, 1.0]])
    with pytest.raises(NotImplementedError, match="SYMMETRIC"):
        MultiSpeciesMetropolis(
            sg,
            species=2,
            q_per_species=3,
            interaction_matrix=M,
            savedisk=False,
        )


def test_unwired_axes(tmp_path):
    sg = make_lattice(tmp_path)
    for move in ("wolff", "sw"):
        with pytest.raises(NotImplementedError, match="cluster moves"):
            MultiSpeciesMetropolis(sg, move=move, savedisk=False)
    with pytest.raises(NotImplementedError, match="conditional resample"):
        MultiSpeciesMetropolis(sg, rule="heatbath", savedisk=False)
    with pytest.raises(NotImplementedError, match="deterministic local"):
        MultiSpeciesMetropolis(sg, upd_mode="gillespie", savedisk=False)


def test_axis_capability_errors(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="species must be >= 1"):
        MultiSpeciesMetropolis(sg, species=0, savedisk=False)
    with pytest.raises(ValueError, match="q_per_species must be >= 2"):
        MultiSpeciesMetropolis(sg, q_per_species=1, savedisk=False)
    with pytest.raises(ValueError, match="length must match"):
        MultiSpeciesMetropolis(
            sg, species=2, q_per_species=(3, 3, 3), savedisk=False
        )
    with pytest.raises(ValueError, match="rule"):
        MultiSpeciesMetropolis(sg, rule="bogus", savedisk=False)
    with pytest.raises(ValueError, match="move"):
        MultiSpeciesMetropolis(sg, move="bogus", savedisk=False)
    with pytest.raises(ValueError, match="asynchronous kinetics"):
        MultiSpeciesMetropolis(
            sg, move="exchange", upd_mode="sync", savedisk=False
        )
    with pytest.raises(ValueError, match="tie_flip_p is metropolis-only"):
        MultiSpeciesMetropolis(
            sg, rule="glauber", tie_flip_p=0.5, savedisk=False
        )
    with pytest.raises(ValueError, match="observable"):
        MultiSpeciesMetropolis(sg, observables=("bogus",), savedisk=False)
    with pytest.raises(ValueError, match="coupling_norm"):
        MultiSpeciesMetropolis(sg, coupling_norm="bogus", savedisk=False)
    with pytest.raises(NotImplementedError, match="field"):
        MultiSpeciesMetropolis(sg, field=np.ones(sg.N), savedisk=False)


def test_pb_capability_guards(tmp_path):
    """The pb solver's supports() gate rejects, at setup, everything the
    compiled kernel cannot faithfully represent (no binary needed)."""
    sg = make_lattice(tmp_path)

    def pb(**kw):
        kw.setdefault("observables", ("energy",))
        kw.setdefault("q_per_species", 3)
        m = MultiSpeciesMetropolis(
            sg,
            species=2,
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
    with pytest.raises(NotImplementedError, match="IDENTITY"):
        pb(interaction_matrix=np.array([[1.0, 0.5], [0.5, 1.0]]))
    with pytest.raises(NotImplementedError, match="UNIFORM"):
        pb(q_per_species=(2, 5))
    with pytest.raises(NotImplementedError, match="energy only"):
        pb(observables=("energy", "magn"))
    with pytest.raises(NotImplementedError, match="snapshots"):
        pb(observables=("energy", "snapshots"))


# =========================================================================
# Native pybind backend (kernel semantics; skipped without the binary)
# =========================================================================
@needs_native
def test_pb_bookkeeping_and_reproducibility(tmp_path):
    sg = make_lattice(tmp_path)
    runs = []
    for _ in range(2):
        m = MultiSpeciesMetropolis(
            sg,
            species=2,
            q_per_species=3,
            T=1.0,
            steps=6,
            seed=SEED,
            savedisk=False,
            runlang="pb",
            observables=("energy",),
        )
        m.run()
        # exact record-then-sweep reconstruction + running-energy invariant
        assert np.isclose(m.ene[0], m.energy_intensive(m.sini))
        assert np.isclose(m._e_running, m.total_energy())
        runs.append((m.s.copy(), np.asarray(m.ene).copy()))
    assert np.array_equal(runs[0][0], runs[1][0])
    assert np.allclose(runs[0][1], runs[1][1])
    other = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=1.0,
        steps=6,
        seed=SEED + 1,
        savedisk=False,
        runlang="pb",
        observables=("energy",),
    )
    other.run()
    assert not np.array_equal(runs[0][0], other.s)


@needs_native
def test_pb_zero_T_is_monotone(tmp_path):
    """The kernel accepts every ΔE <= 0 proposal, so at T = 0 the energy
    trace is nonincreasing (ties flip freely but cost nothing)."""
    sg = make_lattice(tmp_path, pflip=0.1)
    m = MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=0.0,
        steps=20,
        seed=SEED,
        savedisk=False,
        runlang="pb",
        observables=("energy",),
    )
    m.run()
    assert np.all(np.diff(np.asarray(m.ene)) <= 1e-12)


@needs_native
def test_pb_matches_py_equilibrium_energy(tmp_path):
    """Different RNG streams — compare the equilibrium energy statistically
    (D-B6: statistical, never byte-equality)."""

    def tail_energy(runlang, T):
        vals = []
        for sd in (1, 2, 3):
            sg = make_lattice(tmp_path, pflip=0.0)
            m = MultiSpeciesMetropolis(
                sg,
                species=2,
                q_per_species=3,
                T=T,
                steps=150,
                seed=sd,
                savedisk=False,
                runlang=runlang,
                observables=("energy",),
            )
            m.run()
            vals.append(float(np.mean(m.ene[-30:])))
        return float(np.mean(vals))

    for T in (0.5, 3.0):
        assert abs(tail_energy("pb", T) - tail_energy("py", T)) < 0.2
