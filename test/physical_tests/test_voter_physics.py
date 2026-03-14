"""
Test Voter Model physics.

Known behaviors:
- Spin values always in {-1, +1}
- Small unfrustrated graph: consensus (|M| -> 1) over time
- Voter dynamics: magnetization is a martingale (conserved on average)
"""

import pytest
import numpy as np


# ===================================================================
# Spin validity
# ===================================================================


@pytest.mark.physical
def test_voter_spin_values(tmp_path):
    """Voter spins must always be in {-1, +1}."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    lat = Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter = VoterModel(sg=lat, steps=50, runlang="py", seed=42)
    voter.init_voter_dynamics()
    voter.run(tqdm_on=False)

    assert np.all(np.isin(voter.s, [-1, 1])), (
        f"Voter spins should be +/-1, got unique={np.unique(voter.s)}"
    )


# ===================================================================
# Consensus on small graphs
# ===================================================================


@pytest.mark.physical
def test_voter_consensus_complete_graph(tmp_path):
    """Small complete unfrustrated graph should reach consensus."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    from lrgsglib.statsys import VoterModel

    g = FullyConnectedNX(
        N=8, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter = VoterModel(sg=g, steps=500, runlang="py", seed=42)
    voter.init_voter_dynamics()
    voter.run(tqdm_on=False)

    # Complete graph with 8 nodes should reach consensus quickly
    final_mag = np.abs(np.mean(voter.s))
    assert final_mag > 0.8, (
        f"Complete graph should reach near-consensus, got |M|={final_mag:.3f}"
    )


# ===================================================================
# C backend
# ===================================================================


@pytest.mark.physical
@pytest.mark.integration
def test_voter_c_backend_spins(tmp_path):
    """C backend: voter spins must be valid."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    lat = Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter = VoterModel(sg=lat, steps=50, runlang="C0S", seed=42)
    voter.init_voter_dynamics()
    voter.run(tqdm_on=False, verbose=False)

    assert np.all(np.isin(voter.s, [-1, 1]))


# ===================================================================
# Visual: magnetization dynamics
# ===================================================================


@pytest.mark.physical
@pytest.mark.visual
def test_voter_magnetization_visual(tmp_path, output_dir, save_figure, system_size, n_steps):
    """Visual: voter magnetization time series on different graphs."""
    import matplotlib.pyplot as plt
    from lrgsglib.graphs.nx import Lattice2DNX, FullyConnectedNX
    from lrgsglib.statsys import VoterModel

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Small complete graph (should reach consensus fast)
    g_complete = FullyConnectedNX(
        N=12, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter1 = VoterModel(
        sg=g_complete, steps=n_steps, runlang="py", seed=42,
        save_magnetization=True,
    )
    voter1.init_voter_dynamics()
    voter1.run(tqdm_on=False)

    if hasattr(voter1, "magn") and len(voter1.magn) > 0:
        axes[0].plot(voter1.magn)
    axes[0].set_xlabel("Step")
    axes[0].set_ylabel("M(t)")
    axes[0].set_title("Voter on Complete Graph (N=12)")
    axes[0].axhline(0, color="gray", linestyle="--", alpha=0.5)

    # Lattice (slower consensus)
    side = max(system_size, 8)
    lat = Lattice2DNX(
        side1=side, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter2 = VoterModel(
        sg=lat, steps=n_steps, runlang="py", seed=42,
        save_magnetization=True,
    )
    voter2.init_voter_dynamics()
    voter2.run(tqdm_on=False)

    if hasattr(voter2, "magn") and len(voter2.magn) > 0:
        axes[1].plot(voter2.magn)
    axes[1].set_xlabel("Step")
    axes[1].set_ylabel("M(t)")
    axes[1].set_title(f"Voter on 2D Lattice ({side}x{side})")
    axes[1].axhline(0, color="gray", linestyle="--", alpha=0.5)

    fig.suptitle("Voter Model Magnetization Dynamics")
    fig.tight_layout()
    save_figure(fig, output_dir / "voter_magnetization.png")
