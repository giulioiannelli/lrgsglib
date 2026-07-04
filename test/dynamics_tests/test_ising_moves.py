"""Phase-2 Ising Layer-2 moves: wolff / sw / exchange / spectral.

Covers the move axis on IsingMetropolis and the shared engines:

- cluster ΔE bookkeeping: flip_cluster's boundary-bond ΔE must equal the
  total-energy difference (any cluster, any coupling_norm);
- Wolff and Swendsen-Wang reach the exact ground state on the clean
  lattice at low T, and agree with single-site Metropolis in equilibrium
  (statistical check — the detailed-balance smoke test);
- Kawasaki (move='exchange') conserves the magnetization exactly and its
  swap ΔE matches the total-energy difference (including the adjacent-pair
  correction −4 J_uv s_u s_v);
- the spectral move seeds from the best eigenvector, tracks its best-so-far
  consistently, and improves on a frustrated instance;
- axis-interaction capability errors are raised at setup (plan §9).

Citations: Swendsen-Wang (M5), Wolff (M6), Kawasaki (M8) per
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
# Cluster ΔE bookkeeping
# =========================================================================
@pytest.mark.parametrize("coupling_norm", ["raw", "sym"])
def test_flip_cluster_delta_E_matches_energy_difference(
    tmp_path, coupling_norm
):
    sg = make_lattice(tmp_path, pflip=0.15)
    field = np.linspace(-0.5, 0.5, N)
    model = IsingMetropolis(
        sg,
        T=1.0,
        seed=SEED,
        savedisk=False,
        coupling_norm=coupling_norm,
        field=field,
    )
    model.init_dynamics()
    rng = np.random.default_rng(0)
    for size in (1, 5, 20, N):
        nodes = rng.choice(N, size=min(size, N), replace=False)
        E_before = model.total_energy()
        dE = model.flip_cluster(nodes)
        assert dE == pytest.approx(model.total_energy() - E_before, abs=1e-9)


def test_swap_delta_E_matches_energy_difference(tmp_path):
    """Including the adjacent-pair correction −4 J_uv s_u s_v."""
    sg = make_lattice(tmp_path, pflip=0.15)
    model = IsingMetropolis(
        sg,
        T=1.0,
        move="exchange",
        seed=SEED,
        savedisk=False,
        field=np.linspace(-0.5, 0.5, N),
    )
    model.init_dynamics()
    rng = np.random.default_rng(1)
    checked_adjacent = 0
    for _ in range(40):
        u = int(rng.integers(0, N))
        nbrs = model.neighbor_indices(u)
        v = int(nbrs[rng.integers(0, len(nbrs))])
        if model.s[u] == model.s[v]:
            continue
        checked_adjacent += 1
        dE = model.swap_delta_E(u, v)
        swapped = model.s.copy()
        swapped[u], swapped[v] = swapped[v], swapped[u]
        dE_ref = model.total_energy(swapped) - model.total_energy()
        assert dE == pytest.approx(dE_ref, abs=1e-9)
    assert checked_adjacent > 0


# =========================================================================
# Cluster kinetics
# =========================================================================
@pytest.mark.parametrize("move", ["wolff", "sw"])
def test_cluster_moves_reach_ground_state(tmp_path, move):
    """Low T on the clean lattice: the cluster kinetics must find the
    fully-ordered state (e = −2 per site on the periodic square lattice)."""
    sg = make_lattice(tmp_path, pflip=0.0)
    model = IsingMetropolis(
        sg,
        T=0.5,
        move=move,
        steps=20,
        seed=SEED,
        savedisk=False,
        ic="uniform",
    )
    model.run()
    assert model.energy_intensive() == pytest.approx(-2.0)
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-8)
    stats = model.cluster_sizes if move == "wolff" else model.cluster_counts
    assert stats is not None and len(stats) == 20


@pytest.mark.parametrize("move", ["wolff", "sw"])
def test_cluster_agrees_with_single_site_in_equilibrium(tmp_path, move):
    """Same Gibbs state: cluster and single-site kinetics must agree on the
    mean equilibrium energy (loose statistical tolerance)."""
    T = 2.0  # below T_c: ordered, low-variance energies

    def mean_tail_energy(mv, steps, burn):
        sg = make_lattice(tmp_path / mv, pflip=0.0)
        model = IsingMetropolis(
            sg,
            T=T,
            move=mv,
            steps=steps,
            seed=SEED,
            savedisk=False,
            ic="uniform",
        )
        model.run()
        return float(np.mean(model.ene[burn:]))

    e_single = mean_tail_energy("single", 400, 200)
    e_cluster = mean_tail_energy(move, 120, 60)
    assert e_cluster == pytest.approx(e_single, abs=0.1)


# =========================================================================
# Kawasaki exchange
# =========================================================================
def test_exchange_conserves_magnetization(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.1)
    model = IsingMetropolis(
        sg,
        T=1.5,
        move="exchange",
        steps=30,
        seed=SEED,
        savedisk=False,
        ic="uniform",
    )
    model.init_dynamics()
    m0 = model.magnetization()
    model.run()
    assert model.magnetization() == m0
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-8)
    # The magn trace is constant by construction.
    assert np.ptp(np.asarray(model.magn)) == 0.0


# =========================================================================
# Spectral move
# =========================================================================
def test_spectral_move_runs_and_tracks_best(tmp_path):
    sg = make_lattice(tmp_path, pflip=0.15)
    model = IsingMetropolis(
        sg,
        T=0.5,
        move="spectral",
        steps=8,
        spectral_n_modes=6,
        spectral_chunk_size=4,
        spectral_polish=True,
        spectral_polish_sweeps=10,
        seed=SEED,
        savedisk=False,
        ic="uniform",
    )
    model.run()
    assert len(model.ene) == 8
    assert model._e_running == pytest.approx(model.total_energy(), abs=1e-8)
    # Best-so-far bookkeeping: intensive, consistent, no worse than any
    # recorded point of the walk.
    assert model.spectral_best_energy == pytest.approx(
        model.energy_intensive(model.spectral_best_spins), abs=1e-9
    )
    assert model.spectral_best_energy <= min(model.ene) + 1e-9
    assert model.spectral_coeffs.shape == (6,)


def test_spectral_move_seeded_reproducibility(tmp_path):
    def run_once(sub):
        sg = make_lattice(tmp_path / sub, pflip=0.1)
        model = IsingMetropolis(
            sg,
            T=1.0,
            move="spectral",
            steps=5,
            spectral_n_modes=4,
            spectral_polish=False,
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
def test_move_axis_capability_errors(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="rule does not apply"):
        IsingMetropolis(sg, move="wolff", rule="glauber", savedisk=False)
    with pytest.raises(ValueError, match="tie_flip_p does not apply"):
        IsingMetropolis(sg, move="sw", tie_flip_p=0.5, savedisk=False)
    with pytest.raises(ValueError, match="order does not apply"):
        IsingMetropolis(sg, move="spectral", order="typewriter", savedisk=False)
    with pytest.raises(ValueError, match="upd_mode does not apply"):
        IsingMetropolis(sg, move="wolff", upd_mode="sync", savedisk=False)
    with pytest.raises(NotImplementedError, match="ghost-spin"):
        IsingMetropolis(sg, move="wolff", field=np.ones(N), savedisk=False)
    # The exchange move keeps the rule/order/tie axes.
    model = IsingMetropolis(
        sg,
        T=1.0,
        move="exchange",
        rule="glauber",
        order="permutation",
        savedisk=False,
    )
    assert model.move == "exchange"
