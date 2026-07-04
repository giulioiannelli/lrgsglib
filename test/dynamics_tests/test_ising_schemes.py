"""Phase-2 Ising scheme leaves: SA / PT / CEM on the shared spine.

Covers the scheme-level physics and the configure-then-run contract:

- SA: schedule generators (linear/exponential/logarithmic/custom), the
  engine-recompile-per-stage loop cools the system (energy decreases along
  the schedule), running-energy consistency, the spectral (TFCA) field
  enters ΔE and E first-class, fixed step budget (run() rejects steps=);
- PT: ladder generators + validation, per-rung traces, incremental replica
  energies consistent with the final states, the corrected exchange
  criterion min(1, e^{Δβ·ΔE}) — the legacy `_attempt_exchange` accepted
  every swap unconditionally (tautological condition), so the new scheme's
  acceptance must be < 1 somewhere on a frustrated instance;
- CEM: incumbent trace, best-state bookkeeping consistency, the objective
  is the FULL §3b Hamiltonian (field included — exact ground state on the
  unfrustrated lattice with a strong uniform field);
- seeded reproducibility for all three schemes.

Citations: acceptance rules per ``.agents/ref/00-REFERENCES.md`` §MCMC.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import (
    IsingCEM,
    IsingParallelTempering,
    IsingSimulatedAnnealing,
)
from lrgsglib.statsys.IsingDynamics.IsingParallelTempering import (
    generate_temperature_ladder,
)
from lrgsglib.statsys.IsingDynamics.IsingSimulatedAnnealing import (
    generate_temperature_schedule,
)

SIDE = 8
N = SIDE * SIDE
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


# =========================================================================
# Simulated annealing
# =========================================================================
def test_sa_schedule_generators():
    lin = generate_temperature_schedule("linear", 4.0, 1.0, 4, 0.9)
    assert lin == pytest.approx([4.0, 3.0, 2.0, 1.0])
    exp = generate_temperature_schedule("exponential", 2.0, 0.0, 3, 0.5)
    assert exp == pytest.approx([2.0, 1.0, 0.5])  # T_final unused (legacy)
    log = generate_temperature_schedule("logarithmic", 3.0, 0.0, 2, 0.9)
    assert log == pytest.approx(3.0 / np.log([2.0, 3.0]))
    with pytest.raises(ValueError, match="schedule"):
        generate_temperature_schedule("bogus", 1.0, 0.1, 5, 0.9)


def test_sa_cools_and_tracks_energy(tmp_path):
    """Annealing 5 → 0.01 on the clean lattice must lower the energy, and
    the incremental running energy must equal a fresh recompute."""
    sg = make_lattice(tmp_path)
    model = IsingSimulatedAnnealing(
        sg,
        T_init=5.0,
        T_final=0.01,
        schedule="linear",
        n_temperatures=10,
        steps_per_T=5,
        seed=SEED,
        savedisk=False,
        ic="uniform",
    )
    model.run()
    assert model.steps == 10 * 5
    assert len(model.ene) == 10 * 5
    assert len(model.sa_energy) == 10
    assert model.sa_energy[-1] < model.sa_energy[0]
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-8)
    curve = model.sa_curve()
    assert curve["T"] == pytest.approx(model.T_schedule)


def test_sa_custom_schedule_and_fixed_budget(tmp_path):
    sg = make_lattice(tmp_path)
    schedule = np.array([2.0, 1.0, 0.5])
    model = IsingSimulatedAnnealing(
        sg,
        T_schedule=schedule,
        steps_per_T=2,
        seed=SEED,
        savedisk=False,
    )
    assert model.n_temperatures == 3 and model.steps == 6
    with pytest.raises(ValueError, match="schedule"):
        model.run(steps=99)
    with pytest.raises(ValueError):
        IsingSimulatedAnnealing(
            sg, T_schedule=np.array([1.0, -0.5]), savedisk=False
        )


def test_sa_spectral_field_is_first_class(tmp_path):
    """field='spectral' (the legacy TFCA) builds a nonzero field at init
    that enters BOTH ΔE and the tracked energy."""
    sg = make_lattice(tmp_path, pflip=0.1)
    model = IsingSimulatedAnnealing(
        sg,
        T_schedule=np.full(4, 1e-3),
        steps_per_T=5,
        field="spectral",
        spectral_n_modes=6,
        seed=SEED,
        savedisk=False,
        ic="uniform",
    )
    model.init_dynamics()
    assert model.spectral_field is not None
    assert np.any(model.spectral_field != 0.0)
    nd = int(np.argmax(np.abs(model.spectral_field)))
    proposal = model.propose_local(nd)
    dE = model.delta_E(nd, proposal)
    coupling = model.local_field(nd) - float(model.field[nd])
    expected = 2.0 * float(model.s[nd]) * (coupling + float(model.field[nd]))
    assert dE == pytest.approx(expected, abs=1e-9)
    model.run()
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-8)


def test_sa_field_axis_validation(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="field"):
        IsingSimulatedAnnealing(sg, field="bogus", savedisk=False)
    explicit = np.linspace(-1.0, 1.0, N)
    model = IsingSimulatedAnnealing(
        sg,
        field=explicit,
        T_schedule=np.array([1.0]),
        steps_per_T=1,
        savedisk=False,
    )
    assert model.field_mode == "uniform"
    assert model.field == pytest.approx(explicit)


def test_sa_seeded_reproducibility(tmp_path):
    def run_once(sub):
        sg = make_lattice(tmp_path / sub)
        model = IsingSimulatedAnnealing(
            sg,
            T_init=3.0,
            T_final=0.1,
            schedule="linear",
            n_temperatures=5,
            steps_per_T=3,
            seed=SEED,
            savedisk=False,
        )
        model.run()
        return np.array(model.ene), model.s.copy()

    e1, s1 = run_once("a")
    e2, s2 = run_once("b")
    assert e1 == pytest.approx(e2)
    assert np.array_equal(s1, s2)


# =========================================================================
# Parallel tempering
# =========================================================================
def test_pt_ladder_generators():
    geo = generate_temperature_ladder("geometric", 0.5, 4.0, 4)
    assert geo[0] == pytest.approx(0.5) and geo[-1] == pytest.approx(4.0)
    assert np.diff(np.log(geo)) == pytest.approx(np.full(3, np.log(8.0) / 3.0))
    lin = generate_temperature_ladder("linear", 1.0, 3.0, 3)
    assert lin == pytest.approx([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="n_replicas"):
        generate_temperature_ladder("geometric", 0.5, 4.0, 1)
    with pytest.raises(ValueError, match="T_min"):
        generate_temperature_ladder("geometric", 0.0, 4.0, 4)
    with pytest.raises(ValueError, match="ladder"):
        generate_temperature_ladder("bogus", 0.5, 4.0, 4)


def test_pt_custom_ladder_validation(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="ascending"):
        IsingParallelTempering(
            sg, T_ladder=np.array([2.0, 1.0]), savedisk=False
        )
    with pytest.raises(ValueError, match="> 0"):
        IsingParallelTempering(
            sg, T_ladder=np.array([0.0, 1.0]), savedisk=False
        )


def test_pt_runs_and_tracks_replicas(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    model = IsingParallelTempering(
        sg,
        n_replicas=4,
        T_min=0.5,
        T_max=4.0,
        steps_per_exchange=2,
        n_exchanges=12,
        seed=SEED,
        savedisk=False,
        ic="uniform",
    )
    model.run()
    assert model.steps == 12
    assert len(model.ene) == 12
    assert model.pt_energy.shape == (4, 12)
    assert model.pt_magn.shape == (4, 12)
    # Incremental replica energies match a fresh recompute of the finals.
    assert model._e_running == pytest.approx(
        model.total_energy(model.pt_final_states[0]), abs=1e-8
    )
    assert np.array_equal(model.s, model.pt_final_states[0])
    # Even-odd alternation: 12 rounds x alternating 2/1 pairs on 4 rungs.
    assert len(model.pt_exchanges) == 6 * 2 + 6 * 1
    rates = model.exchange_rate()
    assert rates.shape == (3,)
    assert np.all((rates >= 0.0) & (rates <= 1.0))


def test_pt_exchange_criterion_rejects_sometimes(tmp_path):
    """The legacy `_attempt_exchange` condition was a tautology (always
    swapped). The corrected min(1, e^{Δβ·ΔE}) criterion must reject some
    swaps on a frustrated instance with a wide ladder."""
    sg = make_lattice(tmp_path, pflip=0.2)
    model = IsingParallelTempering(
        sg,
        n_replicas=3,
        T_min=0.3,
        T_max=6.0,
        steps_per_exchange=3,
        n_exchanges=40,
        seed=SEED,
        savedisk=False,
        ic="uniform",
    )
    model.run()
    accepted = [acc for (_, _, _, acc) in model.pt_exchanges]
    assert not all(accepted)


def test_pt_seeded_reproducibility(tmp_path):
    def run_once(sub):
        sg = make_lattice(tmp_path / sub, pflip=0.1)
        model = IsingParallelTempering(
            sg,
            n_replicas=3,
            T_min=0.5,
            T_max=3.0,
            steps_per_exchange=2,
            n_exchanges=8,
            seed=SEED,
            savedisk=False,
        )
        model.run()
        return model.pt_energy.copy(), model.s.copy()

    e1, s1 = run_once("a")
    e2, s2 = run_once("b")
    assert e1 == pytest.approx(e2)
    assert np.array_equal(s1, s2)


def test_pt_fixed_budget_and_tie_rules(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingParallelTempering(
        sg, n_replicas=2, n_exchanges=3, savedisk=False
    )
    with pytest.raises(ValueError, match="n_exchanges"):
        model.run(steps=99)
    with pytest.raises(ValueError, match="metropolis-only"):
        IsingParallelTempering(
            sg, rule="glauber", tie_flip_p=0.5, savedisk=False
        )


# =========================================================================
# Cross-entropy method
# =========================================================================
def test_cem_finds_exact_ground_state_with_field(tmp_path):
    """On the unfrustrated periodic square lattice with a strong uniform
    field the ground state is all-up with e = -(2 + h) per site — and the
    CEM objective must be the FULL §3b Hamiltonian (field included)."""
    sg = make_lattice(tmp_path, pflip=0.0)
    h = 2.0
    model = IsingCEM(
        sg,
        field=np.full(N, h),
        n_modes=6,
        pop_size=16,
        n_iter=5,
        restarts=1,
        greedy=True,
        greedy_sweeps=20,
        seed=SEED,
        savedisk=False,
    )
    model.run()
    assert np.all(model.cem_best_spins == 1)
    assert model.cem_best_energy == pytest.approx(-(2.0 + h))


def test_cem_bookkeeping_consistency(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.15)
    model = IsingCEM(
        sg,
        n_modes=8,
        pop_size=12,
        n_iter=4,
        restarts=2,
        greedy=True,
        greedy_sweeps=10,
        seed=SEED,
        savedisk=False,
    )
    model.run()
    assert model.steps == 2 * 4
    assert len(model.ene) == 2 * 4
    # The reported best is the intensive energy of the reported spins,
    # and the trace's minimum can never beat the global best.
    assert model.cem_best_energy == pytest.approx(
        model.energy_intensive(model.cem_best_spins), abs=1e-9
    )
    assert model.cem_restart_energies.shape == (2,)
    assert model.cem_best_energy == pytest.approx(
        model.cem_restart_energies.min()
    )
    assert min(model.ene) >= model.cem_best_energy - 1e-9
    assert np.array_equal(model.s, model.cem_best_spins)


def test_cem_parameter_validation(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="elite_frac"):
        IsingCEM(sg, elite_frac=0.0, savedisk=False)
    with pytest.raises(ValueError, match="pop_size"):
        IsingCEM(sg, pop_size=0, savedisk=False)
    with pytest.raises(ValueError, match="n_iter"):
        IsingCEM(sg, restarts=2, n_iter=0, savedisk=False)
    model = IsingCEM(sg, restarts=1, n_iter=2, savedisk=False)
    with pytest.raises(ValueError, match="n_iter"):
        model.run(steps=5)


def test_cem_seeded_reproducibility(tmp_path):
    def run_once(sub):
        sg = make_lattice(tmp_path / sub, pflip=0.1)
        model = IsingCEM(
            sg,
            n_modes=6,
            pop_size=8,
            n_iter=3,
            restarts=1,
            greedy=False,
            seed=SEED,
            savedisk=False,
        )
        model.run()
        return np.array(model.ene), model.cem_best_spins.copy()

    e1, s1 = run_once("a")
    e2, s2 = run_once("b")
    assert e1 == pytest.approx(e2)
    assert np.array_equal(s1, s2)
