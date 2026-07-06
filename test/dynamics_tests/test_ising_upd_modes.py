"""Phase-2 Ising update modes: sync (sublattice-parallel) + gillespie (BKL).

Covers the upd_mode axis on IsingMetropolis and the shared ThermalEngine:

- sync greedy coloring: the periodic square lattice is bipartite, so the
  coloring must recover the checkerboard — exactly 2 classes that partition
  the sites with no intra-class bond;
- both modes keep the incremental-energy invariant
  (``_e_running == total_energy()`` after every run);
- zero-temperature quenches are monotone in energy, and the rejection-free
  gillespie quench actually freezes (reaches an absorbing state);
- both modes agree with async single-site kinetics on the mean equilibrium
  energy (same Gibbs state, loose statistical tolerance);
- seeded runs reproduce;
- axis-interaction capability errors are raised at setup (plan §9):
  the async 'order' has no meaning under sync/gillespie, and the exchange
  move is asynchronous-only.

Citations: sublattice/checkerboard updates per M4 (Newman & Barkema 1999),
rejection-free continuous-time kinetics per M3 (Bortz-Kalos-Lebowitz 1975),
``.agents/ref/00-REFERENCES.md`` §MCMC.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import IsingMetropolis

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
# Sync: greedy coloring correctness
# =========================================================================
def test_sync_coloring_is_checkerboard_on_square_lattice(tmp_path):
    """Bipartite lattice -> exactly 2 color classes, partitioning the sites
    with no bond inside a class (the frozen-state ΔE validity condition)."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingMetropolis(
        sg, T=1.0, upd_mode="sync", steps=3, seed=SEED, savedisk=False
    )
    model.run()
    classes = model._engine.color_classes
    assert classes is not None and len(classes) == 2
    # Partition: every site exactly once.
    merged = np.sort(np.concatenate(classes))
    assert np.array_equal(merged, np.arange(N))
    # Independence: no site neighbors another site of its own class.
    for cls in classes:
        members = set(int(nd) for nd in cls)
        for nd in members:
            assert not members.intersection(
                int(j) for j in model.neighbor_indices(nd)
            )


# =========================================================================
# Zero-temperature quenches (monotone energy + bookkeeping invariant)
# =========================================================================
@pytest.mark.parametrize("upd_mode", ["sync", "gillespie"])
def test_zero_T_quench_energy_monotone(tmp_path, upd_mode):
    """T=0 with frozen ties: strictly non-increasing energy under both the
    sublattice-parallel and the rejection-free schedule."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingMetropolis(
        sg,
        T=0.0,
        upd_mode=upd_mode,
        tie_flip_p="frozen",
        ic="uniform",
        steps=30,
        seed=3,
        savedisk=False,
    )
    model.run()
    ene = np.asarray(model.ene, dtype=float)
    assert len(ene) == 30
    assert np.all(np.diff(ene) <= 1e-12)
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-8)


def test_gillespie_quench_freezes(tmp_path):
    """Rejection-free kinetics: the T=0 quench must reach an absorbing state
    (no downhill flip left; then Σr = 0 and time just elapses)."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingMetropolis(
        sg,
        T=0.0,
        upd_mode="gillespie",
        tie_flip_p="frozen",
        ic="uniform",
        steps=30,
        seed=3,
        savedisk=False,
    )
    model.run()
    assert model.is_absorbing()
    # The tail of the energy trace is flat once frozen.
    ene = np.asarray(model.ene, dtype=float)
    assert np.ptp(ene[-5:]) == 0.0


# =========================================================================
# Equilibrium agreement with async (same Gibbs state)
# =========================================================================
@pytest.mark.parametrize("upd_mode", ["sync", "gillespie"])
def test_agrees_with_async_in_equilibrium(tmp_path, upd_mode):
    """Sublattice-parallel blocks and BKL events sample the SAME Gibbs
    measure as random-sequential attempts: the mean equilibrium energies
    must agree (loose statistical tolerance)."""
    T = 2.0  # below T_c: ordered, low-variance energies

    def mean_tail_energy(mode, steps, burn):
        sg = make_lattice(tmp_path / mode, pflip=0.0)
        model = IsingMetropolis(
            sg,
            T=T,
            upd_mode=mode,
            steps=steps,
            seed=SEED,
            savedisk=False,
            ic="uniform",
        )
        model.run()
        assert model._e_running == pytest.approx(model.total_energy(), abs=1e-8)
        return float(np.mean(model.ene[burn:]))

    e_async = mean_tail_energy("async", 400, 200)
    e_mode = mean_tail_energy(upd_mode, 400, 200)
    assert e_mode == pytest.approx(e_async, abs=0.1)


# =========================================================================
# Seeded reproducibility
# =========================================================================
@pytest.mark.parametrize("upd_mode", ["sync", "gillespie"])
def test_seeded_reproducibility(tmp_path, upd_mode):
    def run_once(sub):
        sg = make_lattice(tmp_path / sub, pflip=0.1)
        model = IsingMetropolis(
            sg,
            T=2.0,
            upd_mode=upd_mode,
            steps=10,
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
# Axis-interaction capability errors (plan §9: hard, at setup)
# =========================================================================
def test_upd_mode_axis_capability_errors(tmp_path):
    sg = make_lattice(tmp_path)
    # The async site-visit order has no meaning under sync/gillespie.
    with pytest.raises(ValueError, match="order applies to upd_mode='async'"):
        IsingMetropolis(
            sg, T=1.0, upd_mode="sync", order="permutation", savedisk=False
        )
    with pytest.raises(ValueError, match="order applies to upd_mode='async'"):
        IsingMetropolis(
            sg, T=1.0, upd_mode="gillespie", order="typewriter", savedisk=False
        )
    # The exchange move is asynchronous-only (previously silently ignored).
    with pytest.raises(ValueError, match="asynchronous kinetics only"):
        IsingMetropolis(
            sg, T=1.0, move="exchange", upd_mode="sync", savedisk=False
        )
    with pytest.raises(ValueError, match="asynchronous kinetics only"):
        IsingMetropolis(
            sg, T=1.0, move="exchange", upd_mode="gillespie", savedisk=False
        )
    with pytest.raises(ValueError, match="Unknown upd_mode"):
        IsingMetropolis(sg, T=1.0, upd_mode="bogus", savedisk=False)
