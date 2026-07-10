"""Phase-2 native pybind backend (runlang='pb') for the new Ising classes.

The pb solver polymorphs on the instance: IsingMetropolis -> the native
metropolis kernel, IsingSimulatedAnnealing -> the native annealing kernel
(exponential schedule only — the C linear/logarithmic grids differ from the
python generator), IsingParallelTempering -> deliberately refused (the C
exchange criterion is a tautology: every proposed swap is accepted — the
same bug as the legacy python path), IsingCEM -> not wired (the native
objective is pairwise-only, the new class optimizes the full Hamiltonian).

Capability errors are raised at the run() front door BEFORE any output
stream opens and need no built extension; the kernel tests skip when
``_ising_native`` is not built (``make`` in IsingDynamics/ccore/native/).

Citations: ``.agents/ref/00-REFERENCES.md`` §MCMC (M1 Metropolis, M7
Hukushima-Nemoto for the correct PT criterion the refusal protects).
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import (
    IsingCEM,
    IsingMetropolis,
    IsingParallelTempering,
    IsingSimulatedAnnealing,
)
from lrgsglib.statsys.IsingDynamics.IsingSimulatedAnnealing import (
    generate_temperature_schedule,
)

SIDE = 8
N = SIDE * SIDE
SEED = 42


def _native_available() -> bool:
    try:
        from lrgsglib.statsys.IsingDynamics.ccore import (  # noqa: F401
            _ising_native,
        )
    except Exception:
        return False
    return True


needs_native = pytest.mark.skipif(
    not _native_available(),
    reason="_ising_native not built (make in IsingDynamics/ccore/native/)",
)


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
# Capability errors (front door, no built extension needed)
# =========================================================================
def test_pb_capability_errors_metropolis(tmp_path):
    sg = make_lattice(tmp_path)

    def expect_refusal(match, **kw):
        model = IsingMetropolis(sg, runlang="pb", savedisk=False, **kw)
        with pytest.raises(NotImplementedError, match=match):
            model.run()

    expect_refusal("metropolis acceptance only", T=1.0, rule="glauber")
    expect_refusal("order='permutation'", T=1.0, order="permutation")
    expect_refusal("ΔE=0 policy", T=1.0, tie_flip_p=0.5)
    # the T=0 native kernel freezes ties; the standard tie (1.0) mismatches
    expect_refusal("ΔE=0 policy", T=0.0)
    expect_refusal("pairwise-only", T=1.0, field=np.linspace(-1, 1, N))
    expect_refusal(
        "final state only",
        T=1.0,
        observables=("energy", "snapshots"),
    )
    expect_refusal("asynchronous only", T=1.0, upd_mode="sync")
    expect_refusal("move='single' only", T=1.0, move="exchange")


def test_pb_capability_errors_other_schemes(tmp_path):
    sg = make_lattice(tmp_path)
    # PT: deliberate refusal — the C exchange criterion is a tautology.
    model = IsingParallelTempering(
        sg, n_replicas=3, n_exchanges=2, runlang="pb", savedisk=False
    )
    with pytest.raises(NotImplementedError, match="tautology"):
        model.run()
    # CEM: base default (not wired).
    model = IsingCEM(sg, n_iter=2, restarts=1, runlang="pb", savedisk=False)
    with pytest.raises(NotImplementedError, match="not wired"):
        model.run()
    # SA: only the exponential schedule matches the C cooling grid.
    model = IsingSimulatedAnnealing(
        sg,
        schedule="linear",
        T_init=2.0,
        T_final=0.1,
        n_temperatures=4,
        steps_per_T=2,
        runlang="pb",
        savedisk=False,
    )
    with pytest.raises(NotImplementedError, match="exponential"):
        model.run()
    # SA: the spectral (TFCA) field would be dropped from the native trace.
    model = IsingSimulatedAnnealing(
        sg,
        field="spectral",
        n_temperatures=4,
        steps_per_T=2,
        runlang="pb",
        savedisk=False,
    )
    with pytest.raises(NotImplementedError, match="spectral"):
        model.run()


# =========================================================================
# Native Metropolis kernel
# =========================================================================
@needs_native
def test_pb_met_runs_and_tracks_energy(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    model = IsingMetropolis(
        sg, T=1.5, steps=10, seed=SEED, runlang="pb", savedisk=False
    )
    model.run()
    ene = np.asarray(model.ene, dtype=float)
    magn = np.asarray(model.magn, dtype=float)
    assert ene.shape == (10,) and magn.shape == (10,)
    # trace is recorded before each sweep: ene[0] is the initial state
    assert ene[0] == pytest.approx(model.energy_intensive(model.sini))
    # the incremental-energy invariant holds after the native run too
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-9)


@needs_native
@pytest.mark.parametrize("order", ["random", "typewriter"])
def test_pb_met_seeded_reproducibility(tmp_path, order):
    def run_once(sub):
        sg = make_lattice(tmp_path / sub, pflip=0.1)
        model = IsingMetropolis(
            sg,
            T=2.0,
            steps=10,
            order=order,
            seed=SEED,
            runlang="pb",
            savedisk=False,
        )
        model.run()
        return np.array(model.ene), model.s.copy()

    e1, s1 = run_once("a")
    e2, s2 = run_once("b")
    assert e1 == pytest.approx(e2)
    assert np.array_equal(s1, s2)


@needs_native
def test_pb_met_agrees_with_py_in_equilibrium(tmp_path):
    """Same Gibbs state as the python engine (statistical check; the RNG
    streams differ — SFMT vs numpy — so only the means can be compared).
    Averaged over a few seeds and judged against the seed-to-seed
    scatter: a single-run comparison at a fixed abs tolerance sits on
    the tail of the trajectory-to-trajectory distribution and flakes."""
    T = 2.0
    seeds = (SEED, SEED + 1, SEED + 2)

    def tail_energies(runlang):
        vals = []
        for seed in seeds:
            sg = make_lattice(tmp_path / f"{runlang}{seed}", pflip=0.0)
            model = IsingMetropolis(
                sg,
                T=T,
                steps=400,
                seed=seed,
                runlang=runlang,
                savedisk=False,
                ic="uniform",
            )
            model.run()
            vals.append(float(np.mean(model.ene[200:])))
        return np.asarray(vals)

    e_pb, e_py = tail_energies("pb"), tail_energies("py")
    scatter = float(np.hypot(e_pb.std(), e_py.std()) / np.sqrt(len(seeds)))
    assert abs(e_pb.mean() - e_py.mean()) <= max(0.1, 4.0 * scatter)


@needs_native
def test_pb_met_zero_T_quench_monotone(tmp_path):
    """The native T=0 kernel is the frozen-tie quench (only downhill
    flips), so tie_flip_p='frozen' is the matching python-side policy."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingMetropolis(
        sg,
        T=0.0,
        tie_flip_p="frozen",
        ic="uniform",
        steps=30,
        seed=3,
        runlang="pb",
        savedisk=False,
    )
    model.run()
    ene = np.asarray(model.ene, dtype=float)
    assert np.all(np.diff(ene) <= 1e-12)


# =========================================================================
# Native annealing kernel
# =========================================================================
@needs_native
def test_pb_sa_schedule_matches_python_generator(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    model = IsingSimulatedAnnealing(
        sg,
        T_init=3.0,
        cooling_rate=0.8,
        n_temperatures=12,
        steps_per_T=5,
        seed=SEED,
        runlang="pb",
        savedisk=False,
    )
    model.run()
    expected = generate_temperature_schedule("exponential", 3.0, 0.0, 12, 0.8)
    assert model.sa_temps == pytest.approx(expected)
    # full budget ran (no C-side early exit), sweep-resolution traces
    assert len(model.ene) == 12 * 5
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-9)
    # end-of-stage curve bookkeeping
    assert model.sa_energy.shape == (12,)
    assert model.sa_magn.shape == (12,)
    assert model.sa_energy[-1] == pytest.approx(model.energy_intensive())
    assert model.sa_magn[-1] == pytest.approx(model.magnetization())


@needs_native
def test_pb_sa_reaches_ground_state(tmp_path):
    """Clean lattice, slow-enough cooling: the native anneal must find the
    fully-ordered state (e = -2 per site on the periodic square lattice)."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingSimulatedAnnealing(
        sg,
        T_init=3.0,
        cooling_rate=0.8,
        n_temperatures=25,
        steps_per_T=20,
        seed=SEED,
        runlang="pb",
        savedisk=False,
        ic="uniform",
    )
    model.run()
    assert model.energy_intensive() == pytest.approx(-2.0)
