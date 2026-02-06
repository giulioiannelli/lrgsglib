"""
Test Contact Process physics in known regimes.

Key behaviors:
- Absorbing state (all-inactive) must remain absorbing
- Subcritical (gamma << gamma_c): extinction (rho -> 0)
- Supercritical (gamma >> gamma_c): survival (rho > 0)
- 2D without negative links: bistable (no power-law)
- 3D: power-law decay possible even without negative links
- Near critical with negative edges: fractal activity cascades
"""

import pytest
import numpy as np


# ===================================================================
# Absorbing state
# ===================================================================


@pytest.mark.physical
def test_cp_density_bounds(tmp_path):
    """Contact process density must be bounded: 0 <= rho <= 1."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    lat = Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    cp = ContactProcessEI(
        sg=lat, gamma=1.5, steps=50, runlang="py", seed=42,
    )
    cp.init_contact_dynamics()
    cp.run(tqdm_on=False)

    # Final state values must be in {0, 1} (binary) or {-1, 1} (bipolar)
    unique_vals = np.unique(cp.s)
    valid_bipolar = np.all(np.isin(unique_vals, [-1, 1]))
    valid_binary = np.all(np.isin(unique_vals, [0, 1]))
    assert valid_bipolar or valid_binary, (
        f"CP state values should be {{0,1}} or {{-1,1}}, got {unique_vals}"
    )


# ===================================================================
# Phase behavior: extreme regimes
# ===================================================================


@pytest.mark.physical
def test_cp_supercritical_survival(tmp_path, n_steps):
    """Supercritical gamma >> gamma_c: density should remain > 0."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    cp = ContactProcessEI(
        sg=lat, gamma=5.0, steps=n_steps, runlang="py", seed=42,
    )
    cp.init_contact_dynamics()
    cp.run(tqdm_on=False)

    # Supercritical: most nodes should be active
    final_density = np.mean(cp.s == 1) if 0 in cp.s else np.mean(cp.s > 0)
    assert final_density > 0.3, (
        f"Supercritical (gamma=5.0) should survive, got rho={final_density:.3f}"
    )


# ===================================================================
# C backend phase tests
# ===================================================================


@pytest.mark.physical
@pytest.mark.integration
def test_cp_c_supercritical(tmp_path, n_steps):
    """C backend: supercritical gamma must maintain activity."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    cp = ContactProcessEI(
        sg=lat, gamma=5.0, steps=n_steps, runlang="C1c", seed=42,
    )
    cp.init_contact_dynamics()
    cp.run(tqdm_on=False, verbose=False)

    final_density = np.mean(cp.s == 1) if 0 in cp.s else np.mean(cp.s > 0)
    assert final_density > 0.3, (
        f"C backend supercritical should survive, got rho={final_density:.3f}"
    )


# ===================================================================
# Visual: phase diagram and density time series
# ===================================================================


@pytest.mark.physical
@pytest.mark.visual
def test_cp_density_vs_gamma(tmp_path, output_dir, save_figure, system_size):
    """Visual: density vs gamma showing subcritical/supercritical regimes."""
    import matplotlib.pyplot as plt
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    gammas = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
    side = max(system_size, 8)

    lat = Lattice2DNX(
        side1=side, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )

    final_densities = []
    for gamma in gammas:
        cp = ContactProcessEI(
            sg=lat, gamma=gamma, steps=200, runlang="py", seed=42,
        )
        cp.init_contact_dynamics()
        cp.run(tqdm_on=False)
        rho = np.mean(cp.s == 1) if 0 in cp.s else np.mean(cp.s > 0)
        final_densities.append(rho)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(gammas, final_densities, "o-", markersize=8)
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$\rho_\infty$")
    ax.set_title(f"Contact Process Phase Diagram (2D sqr, N={side}x{side})")
    ax.axhline(0, color="gray", linestyle="--", alpha=0.5)

    fig.tight_layout()
    save_figure(fig, output_dir / "cp_phase_diagram.png")


@pytest.mark.physical
@pytest.mark.visual
def test_cp_density_time_series(tmp_path, output_dir, save_figure, system_size):
    """Visual: density rho(t) for several gamma values (manual tracking)."""
    import matplotlib.pyplot as plt
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    side = max(system_size, 8)
    n_snapshots = 20

    fig, ax = plt.subplots(figsize=(8, 5))

    for gamma in [0.5, 1.5, 3.0, 5.0]:
        densities = []
        for step_chunk in range(n_snapshots):
            lat = Lattice2DNX(
                side1=side, geo="sqr", pbc=True, pflip=0.0,
                seed=42, path_data=tmp_path, path_plot=tmp_path,
                init_nw_dict=False,
            )
            cp = ContactProcessEI(
                sg=lat, gamma=gamma, steps=(step_chunk + 1) * 15,
                runlang="py", seed=42,
            )
            cp.init_contact_dynamics()
            cp.run(tqdm_on=False)
            rho = np.mean(cp.s == 1) if 0 in cp.s else np.mean(cp.s > 0)
            densities.append(rho)
        ax.plot(densities, label=f"$\\gamma={gamma}$")

    ax.set_xlabel("Snapshot")
    ax.set_ylabel(r"$\rho$")
    ax.set_title("Contact Process Density Time Series")
    ax.legend()

    fig.tight_layout()
    save_figure(fig, output_dir / "cp_density_timeseries.png")


# ===================================================================
# 2D vs 3D comparison
# ===================================================================


@pytest.mark.physical
@pytest.mark.visual
@pytest.mark.slow
def test_cp_2d_vs_3d(tmp_path, output_dir, save_figure):
    """Visual: 2D (bistable) vs 3D (power-law possible) comparison."""
    import matplotlib.pyplot as plt
    from lrgsglib.graphs.nx import Lattice2DNX, Lattice3DNX
    from lrgsglib.statsys import ContactProcessEI

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 2D
    lat2d = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    # 2D: run multiple short sims to build density curve
    densities_2d = []
    for t in range(1, 21):
        cp2d = ContactProcessEI(
            sg=lat2d, gamma=1.7, steps=t * 25, runlang="py", seed=42,
        )
        cp2d.init_contact_dynamics()
        cp2d.run(tqdm_on=False)
        rho = np.mean(cp2d.s == 1) if 0 in cp2d.s else np.mean(cp2d.s > 0)
        densities_2d.append(rho)

    axes[0].plot(densities_2d)
    axes[0].set_xlabel("Snapshot")
    axes[0].set_ylabel(r"$\rho$")
    axes[0].set_title("2D (pflip=0, bistable)")

    # 3D
    lat3d = Lattice3DNX(
        dim=6, geo="sc", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    densities_3d = []
    for t in range(1, 21):
        cp3d = ContactProcessEI(
            sg=lat3d, gamma=1.7, steps=t * 25, runlang="py", seed=42,
        )
        cp3d.init_contact_dynamics()
        cp3d.run(tqdm_on=False)
        rho = np.mean(cp3d.s == 1) if 0 in cp3d.s else np.mean(cp3d.s > 0)
        densities_3d.append(rho)

    axes[1].plot(densities_3d)
    axes[1].set_xlabel("Snapshot")
    axes[1].set_ylabel(r"$\rho$")
    axes[1].set_title("3D SC (pflip=0)")

    fig.suptitle(r"Contact Process: 2D vs 3D at $\gamma=1.7$")
    fig.tight_layout()
    save_figure(fig, output_dir / "cp_2d_vs_3d.png")


# ===================================================================
# Frustrated CP: signed edges
# ===================================================================


@pytest.mark.physical
@pytest.mark.visual
def test_cp_frustrated_dynamics(tmp_path, output_dir, save_figure, system_size, n_steps):
    """Visual: CP on frustrated lattice with negative edges."""
    import matplotlib.pyplot as plt
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    side = max(system_size, 8)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, pflip, label in zip(axes, [0.0, 0.3], ["pflip=0", "pflip=0.3"]):
        densities = []
        for t in range(1, 16):
            lat = Lattice2DNX(
                side1=side, geo="sqr", pbc=True, pflip=pflip,
                seed=42, path_data=tmp_path, path_plot=tmp_path,
                init_nw_dict=False,
            )
            if pflip > 0:
                lat.flip_random_fract_edges()

            cp = ContactProcessEI(
                sg=lat, gamma=2.0, steps=t * (n_steps // 15 + 1),
                runlang="py", seed=42,
            )
            cp.init_contact_dynamics()
            cp.run(tqdm_on=False)
            rho = np.mean(cp.s == 1) if 0 in cp.s else np.mean(cp.s > 0)
            densities.append(rho)

        ax.plot(densities)
        ax.set_xlabel("Snapshot")
        ax.set_ylabel(r"$\rho$")
        ax.set_title(f"CP {label}, $\\gamma=2.0$")

    fig.suptitle("Effect of Frustration on Contact Process")
    fig.tight_layout()
    save_figure(fig, output_dir / "cp_frustrated_comparison.png")
