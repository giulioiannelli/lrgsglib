"""
Test Ising model physics in known regimes.

Uses extreme parameter values where behavior is deterministic:
- T -> 0: ground state, |M| -> 1 (ferromagnet)
- T -> inf: disorder, |M| -> 0
- T = T_c ~ 2.269: critical (visual inspection)
- pflip > 0.5: frustrated/glassy (visual inspection)
"""

import pytest
import numpy as np


# ===================================================================
# Extreme temperature regimes (deterministic)
# ===================================================================


@pytest.mark.physical
def test_ising_ground_state(tmp_path, n_steps):
    """T -> 0: ferromagnet must order (|M| -> 1)."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(
        sg=lat, T=0.1, steps=n_steps, runlang="py", seed=42,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)

    final_mag = np.abs(np.mean(ising.s))
    assert final_mag > 0.9, f"T->0 should give |M|->1, got {final_mag:.3f}"


@pytest.mark.physical
def test_ising_high_temperature(tmp_path, n_steps):
    """T -> inf: system must disorder (|M| -> 0)."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(
        sg=lat, T=100.0, steps=n_steps, runlang="py", seed=42,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)

    final_mag = np.abs(np.mean(ising.s))
    assert final_mag < 0.3, f"T->inf should give |M|->0, got {final_mag:.3f}"


# ===================================================================
# Spin and energy bounds
# ===================================================================


@pytest.mark.physical
def test_ising_valid_spins(tmp_path):
    """Spin values must always be in {-1, +1}."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(
        sg=lat, T=2.0, steps=20, runlang="py", seed=42,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)

    assert np.all(np.isin(ising.s, [-1, 1])), "All spins must be +/-1"


@pytest.mark.physical
def test_ising_energy_bounds(tmp_path):
    """Ising energy must be bounded: -Ne <= E <= Ne."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(
        sg=lat, T=2.0, steps=20, runlang="py", seed=42,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)

    E = ising.calc_full_energy()
    Ne = lat.Ne
    assert -Ne <= E <= Ne, f"Energy {E} outside bounds [-{Ne}, {Ne}]"


# ===================================================================
# Frustration effects
# ===================================================================


@pytest.mark.physical
def test_frustration_increases_energy(tmp_path, n_steps):
    """Frustrated lattice has higher ground state energy than unfrustrated."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    # Unfrustrated
    lat_uf = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising_uf = IsingDynamics(
        sg=lat_uf, T=0.1, steps=n_steps, runlang="py", seed=42,
    )
    ising_uf.init_ising_dynamics()
    ising_uf.run(tqdm_on=False)
    E_uf = ising_uf.calc_full_energy()

    # Frustrated
    lat_f = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.5,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    lat_f.flip_random_fract_edges()
    ising_f = IsingDynamics(
        sg=lat_f, T=0.1, steps=n_steps, runlang="py", seed=42,
    )
    ising_f.init_ising_dynamics()
    ising_f.run(tqdm_on=False)
    E_f = ising_f.calc_full_energy()

    assert E_f > E_uf, (
        f"Frustrated E={E_f} should be higher than unfrustrated E={E_uf}"
    )


@pytest.mark.physical
def test_frustrated_lattice_no_full_ordering(tmp_path, n_steps):
    """Frustrated lattice (pflip=0.6) cannot fully order even at low T."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.6,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    lat.flip_random_fract_edges()

    ising = IsingDynamics(
        sg=lat, T=0.5, steps=n_steps, runlang="py", seed=42,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)

    final_mag = np.abs(np.mean(ising.s))
    assert final_mag < 0.8, (
        f"Frustrated system should not fully order, got |M|={final_mag:.3f}"
    )


# ===================================================================
# Seed reproducibility
# ===================================================================


@pytest.mark.physical
def test_ising_seed_reproducibility(tmp_path):
    """Same seed must produce identical results."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    results = []
    for _ in range(2):
        lat = Lattice2DNX(
            side1=4, geo="sqr", pbc=True, pflip=0.0,
            seed=42, path_data=tmp_path, path_plot=tmp_path,
            init_nw_dict=False,
        )
        ising = IsingDynamics(
            sg=lat, T=2.0, steps=10, runlang="py", seed=123,
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)
        results.append(ising.s.copy())

    np.testing.assert_array_equal(results[0], results[1])


# ===================================================================
# C backend physics (same invariants)
# ===================================================================


@pytest.mark.physical
@pytest.mark.integration
def test_ising_c_ground_state(tmp_path, n_steps):
    """C backend: T -> 0 must give ordered state."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    # C backend needs more sweeps than Python to order
    c_steps = max(n_steps, 500)
    ising = IsingDynamics(
        sg=lat, T=0.1, steps=c_steps, runlang="C0E", seed=42,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False, verbose=False)

    final_mag = np.abs(np.mean(ising.s))
    assert final_mag > 0.9, f"C backend T->0 should give |M|->1, got {final_mag:.3f}"


@pytest.mark.physical
@pytest.mark.integration
def test_ising_c_high_temperature(tmp_path, n_steps):
    """C backend: T -> inf must give disordered state."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(
        sg=lat, T=100.0, steps=n_steps, runlang="C0E", seed=42,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False, verbose=False)

    final_mag = np.abs(np.mean(ising.s))
    assert final_mag < 0.3, f"C backend T->inf should give |M|->0, got {final_mag:.3f}"


# ===================================================================
# Visual: Critical temperature and frustrated dynamics
# ===================================================================


@pytest.mark.physical
@pytest.mark.visual
@pytest.mark.slow
def test_ising_critical_temperature_visual(tmp_path, output_dir, save_figure, system_size):
    """Visual: Ising at T_c shows domains at all scales."""
    import matplotlib.pyplot as plt
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    T_c = 2.0 / np.log(1 + np.sqrt(2))  # Exact: ~2.269
    side = max(system_size, 16)  # Need at least 16 for visual

    lat = Lattice2DNX(
        side1=side, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(
        sg=lat, T=T_c, steps=1000, runlang="py", seed=42,
        save_magnetization=True,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Spin configuration
    axes[0].imshow(ising.s.reshape(side, side), cmap="coolwarm", vmin=-1, vmax=1)
    axes[0].set_title(f"Ising at $T_c = {T_c:.3f}$, N={side}x{side}")

    # Magnetization time series
    axes[1].plot(ising.magn)
    axes[1].set_xlabel("Step")
    axes[1].set_ylabel("M(t)")
    axes[1].set_title("Magnetization at $T_c$")
    axes[1].axhline(0, color="gray", linestyle="--", alpha=0.5)

    fig.tight_layout()
    save_figure(fig, output_dir / "ising_critical_T.png")


@pytest.mark.physical
@pytest.mark.visual
def test_ising_frustrated_dynamics_visual(tmp_path, output_dir, save_figure, system_size, n_steps):
    """Visual: Frustrated lattice shows glassy dynamics."""
    import matplotlib.pyplot as plt
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    side = max(system_size, 8)
    lat = Lattice2DNX(
        side1=side, geo="sqr", pbc=True, pflip=0.6,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    lat.flip_random_fract_edges()

    ising = IsingDynamics(
        sg=lat, T=0.5, steps=n_steps, runlang="py", seed=42,
        save_magnetization=True,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].imshow(ising.s.reshape(side, side), cmap="coolwarm", vmin=-1, vmax=1)
    axes[0].set_title(f"Frustrated Ising (pflip=0.6, T=0.5)")

    axes[1].plot(ising.magn)
    axes[1].set_xlabel("Step")
    axes[1].set_ylabel("M(t)")
    axes[1].set_title("Magnetization (glassy dynamics)")
    axes[1].axhline(0, color="gray", linestyle="--", alpha=0.5)

    fig.tight_layout()
    save_figure(fig, output_dir / "ising_frustrated_dynamics.png")
