"""Phase-1 Ising vertical slice: IsingMetropolis on the shared spine.

Covers the §3b physics fixes and the configure-then-run contract:

- ΔE from the substrate matches the total-energy difference for EVERY
  coupling_norm (the matched-J fix — the legacy python kernel violated this);
- exact ground-state energy on the unfrustrated lattice;
- zero-temperature quenches are monotone in energy (metropolis and glauber);
- tie presets / capability errors are raised at setup, never mid-run;
- the run() front door records + persists through the ObservableSet;
- the legacy IsingDynamics python path is untouched.

Citations: acceptance rules per ``.agents/ref/00-REFERENCES.md`` §MCMC
(M1 Metropolis 1953, M2 Glauber 1963).
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import IsingDynamics, IsingMetropolis
from lrgsglib.statsys.IsingDynamics.defaults import ISING_COUPLING_NORMS

SIDE = 8
N = SIDE * SIDE
STEPS = 30
SEED = 42


def make_lattice(tmp_path, pflip=0.0):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=7,
    )


@pytest.mark.parametrize("coupling_norm", ISING_COUPLING_NORMS)
@pytest.mark.parametrize("pflip", [0.0, 0.15])
def test_delta_E_matches_energy_difference(tmp_path, coupling_norm, pflip):
    """The single symmetric J feeds both ΔE and E (plan §3b)."""
    sg = make_lattice(tmp_path, pflip=pflip)
    model = IsingMetropolis(
        sg,
        T=1.0,
        seed=SEED,
        savedisk=False,
        coupling_norm=coupling_norm,
    )
    model.init_dynamics()
    rng = np.random.default_rng(0)
    for nd in rng.integers(0, model.N, size=20):
        nd = int(nd)
        proposal = model.propose_local(nd)
        dE = model.delta_E(nd, proposal)
        flipped = model.s.copy()
        flipped[nd] = proposal
        dE_ref = model.total_energy(flipped) - model.total_energy()
        assert dE == pytest.approx(dE_ref, abs=1e-9)


def test_field_enters_delta_E_and_energy(tmp_path):
    """h is first-class in BOTH ΔE and E (legacy python dropped it)."""
    sg = make_lattice(tmp_path)
    field = np.linspace(-1.0, 1.0, N)
    model = IsingMetropolis(sg, T=1.0, seed=SEED, savedisk=False, field=field)
    model.init_dynamics()
    nd = 5
    proposal = model.propose_local(nd)
    dE = model.delta_E(nd, proposal)
    flipped = model.s.copy()
    flipped[nd] = proposal
    assert dE == pytest.approx(
        model.total_energy(flipped) - model.total_energy(), abs=1e-9
    )
    # decomposition: ΔE = 2 s_nd (coupling_nd + h_nd) — the field term is there
    coupling = model.local_field(nd) - field[nd]
    assert dE == pytest.approx(
        2.0 * float(model.s[nd]) * (coupling + field[nd]), abs=1e-9
    )
    assert field[nd] != 0.0


def test_ground_state_energy_exact(tmp_path):
    """All-up on the unfrustrated lattice: E = -Ne (raw J, h = 0), frozen."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingMetropolis(
        sg,
        T=0.0,
        ic="homogeneous",
        seed=SEED,
        savedisk=False,
    )
    model.init_dynamics()
    assert model.total_energy() == pytest.approx(-sg.Ne, abs=1e-9)
    assert model.energy_intensive() == pytest.approx(-sg.Ne / N, abs=1e-12)
    assert model.is_absorbing()


@pytest.mark.parametrize("rule", ["metropolis", "glauber"])
def test_zero_T_quench_energy_monotone(tmp_path, rule):
    """T=0 dynamics never accepts an uphill move (M1/M2 zero-T limits)."""
    sg = make_lattice(tmp_path, pflip=0.0)
    tie = {"rule": rule}
    if rule == "metropolis":
        tie["tie_flip_p"] = "frozen"  # pure downhill: strictly non-increasing
    model = IsingMetropolis(
        sg,
        T=0.0,
        ic="uniform",
        seed=3,
        savedisk=False,
        steps=STEPS,
        **tie,
    )
    model.run()
    ene = np.asarray(model.ene, dtype=float)
    assert len(ene) == STEPS
    assert np.all(np.diff(ene) <= 1e-12)


def test_run_records_and_persists(tmp_path):
    """run() front door: records energy/magn/snapshots, files under dynpath,
    lazy views reload, and no dirs exist before the run (lazy contract)."""
    sg = make_lattice(tmp_path, pflip=0.1)
    model = IsingMetropolis(
        sg,
        T=2.0,
        seed=5,
        savedisk=True,
        steps=10,
        observables=("energy", "magn", "snapshots"),
    )
    assert [p for p in tmp_path.rglob("*") if p.is_dir()] == []
    model.run()
    ene = np.asarray(model.ene, dtype=float)
    magn = np.asarray(model.magn, dtype=float)
    assert ene.shape == (10,) and magn.shape == (10,)
    # incremental energy bookkeeping agrees with a fresh full recompute
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-6)
    traj = np.asarray(model.s_t)
    assert traj.shape == (10, N)
    files = sorted(p.name for p in model.dynpath.iterdir())
    assert any(f.startswith("ene") for f in files)
    assert any(f.startswith("m_") for f in files)
    assert any(f.startswith("sout") for f in files)


def test_seeded_reproducibility(tmp_path):
    finals = []
    for sub in ("a", "b"):
        d = tmp_path / sub
        d.mkdir()
        sg = make_lattice(d, pflip=0.1)
        model = IsingMetropolis(sg, T=2.0, seed=SEED, savedisk=False, steps=10)
        model.run()
        finals.append(model.s.copy())
    assert np.array_equal(finals[0], finals[1])


def test_capability_errors_at_setup(tmp_path):
    sg = make_lattice(tmp_path)
    # upd_mode is wired since Phase 2; the async-only 'order' axis is its
    # capability surface now (full coverage in test_ising_upd_modes.py).
    with pytest.raises(ValueError):
        IsingMetropolis(sg, T=1.0, upd_mode="sync", order="permutation")
    # move='wolff' is wired since Phase 2; its axis interactions are the
    # capability surface now (full coverage in test_ising_moves.py).
    with pytest.raises(ValueError):
        IsingMetropolis(sg, T=1.0, move="wolff", rule="glauber")
    with pytest.raises(ValueError):
        IsingMetropolis(sg, T=1.0, rule="glauber", tie_flip_p=0.5)
    with pytest.raises(ValueError):
        IsingMetropolis(sg, T=1.0, rule="nope")
    with pytest.raises(ValueError):
        IsingMetropolis(sg, T=1.0, tie_flip_p="unknown_preset")
    with pytest.raises(ValueError):
        IsingMetropolis(sg, T=1.0, observables=("energy", "bogus"))
    # native backends: what the compiled kernel cannot represent is a hard
    # error at the front door, before any output opens (full pb coverage in
    # test_ising_pb.py; runlang='cu' stays unwired for the new classes).
    model = IsingMetropolis(sg, T=1.0, rule="glauber", runlang="pb")
    with pytest.raises(NotImplementedError):
        model.run()
    model = IsingMetropolis(sg, T=1.0, runlang="cu")
    with pytest.raises(NotImplementedError):
        model.run()


def test_tie_frozen_is_absorbing_all_up(tmp_path):
    """tie_flip_p='frozen': the all-up unfrustrated state cannot move at T=0."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingMetropolis(
        sg,
        T=0.0,
        ic="homogeneous",
        tie_flip_p="frozen",
        seed=1,
        savedisk=False,
        steps=5,
    )
    model.run()
    assert np.all(model.s == 1)
    assert np.asarray(model.ene, dtype=float) == pytest.approx(-sg.Ne / N)


def test_legacy_ising_python_untouched(tmp_path):
    """Strangler contract: the god-object's python path behaves as before."""
    sg = make_lattice(tmp_path, pflip=0.1)
    legacy = IsingDynamics(sg, T=1.0, runlang="py", eqSTEP=3)
    legacy.init_ising_dynamics()
    legacy.run(tqdm_on=False)
    assert len(legacy.ene) == 3
    assert len(legacy.magn) == 3
