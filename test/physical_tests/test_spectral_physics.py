"""
Test spectral invariants: Laplacian properties, eigenvalue bounds,
backend consistency, entropy/specific-heat observables, DOS, and IPR.

The signed Laplacian is constructed to have non-negative eigenvalues.
"""

import pytest
import numpy as np


# ===================================================================
# Laplacian invariants
# ===================================================================


@pytest.mark.physical
def test_eigenvalues_are_real(small_lattice_sqr):
    """All Laplacian eigenvalues must be real."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    _, eigvals = get_graph_lspectrum(small_lattice_sqr.gr["G"], library="numpy")
    assert np.all(np.isreal(eigvals)), "Eigenvalues must be real"


@pytest.mark.physical
def test_eigenvalues_non_negative(small_lattice_sqr):
    """All Laplacian eigenvalues must be >= 0."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    _, eigvals = get_graph_lspectrum(small_lattice_sqr.gr["G"], library="numpy")
    assert np.all(eigvals >= -1e-10), (
        f"Eigenvalues must be non-negative, min={eigvals.min()}"
    )


@pytest.mark.physical
def test_zero_eigenvalue_connected(small_lattice_sqr):
    """Connected graph: smallest eigenvalue must be ~0."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    _, eigvals = get_graph_lspectrum(small_lattice_sqr.gr["G"], library="numpy")
    eigvals_sorted = np.sort(eigvals)
    assert np.isclose(eigvals_sorted[0], 0.0, atol=1e-10), (
        f"Smallest eigenvalue should be ~0, got {eigvals_sorted[0]}"
    )


@pytest.mark.physical
def test_spectral_gap_connected(small_lattice_sqr):
    """Connected graph: second eigenvalue (spectral gap) must be > 0."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    _, eigvals = get_graph_lspectrum(small_lattice_sqr.gr["G"], library="numpy")
    eigvals_sorted = np.sort(eigvals)
    assert eigvals_sorted[1] > 1e-10, (
        f"Spectral gap should be > 0, got lambda_2={eigvals_sorted[1]}"
    )


@pytest.mark.physical
def test_trace_equals_twice_edges(small_lattice_sqr):
    """Tr(L) = 2 * number of edges."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    L, eigvals = get_graph_lspectrum(small_lattice_sqr.gr["G"], library="numpy")
    trace_L = np.trace(L)
    assert np.isclose(trace_L, 2 * small_lattice_sqr.Ne, rtol=1e-10), (
        f"Tr(L)={trace_L} != 2*Ne={2 * small_lattice_sqr.Ne}"
    )


@pytest.mark.physical
def test_eigenvalue_sum_equals_trace(small_lattice_sqr):
    """Sum of eigenvalues must equal trace of Laplacian."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    L, eigvals = get_graph_lspectrum(small_lattice_sqr.gr["G"], library="numpy")
    assert np.isclose(np.sum(eigvals), np.trace(L), rtol=1e-10)


# ===================================================================
# Backend consistency
# ===================================================================


@pytest.mark.physical
def test_numpy_scipy_consistency(small_lattice_sqr):
    """numpy and scipy backends must produce the same eigenvalues."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    G = small_lattice_sqr.gr["G"]
    _, eigvals_np = get_graph_lspectrum(G, library="numpy")
    _, eigvals_sp = get_graph_lspectrum(G, library="scipy")
    np.testing.assert_allclose(
        np.sort(eigvals_np), np.sort(eigvals_sp), rtol=1e-8,
        err_msg="numpy and scipy eigenvalues differ"
    )


# ===================================================================
# Entropy and specific heat
# ===================================================================


@pytest.mark.physical
def test_entropy_bounds(small_mcg):
    """Entropy must be bounded: 0 <= S <= log(N)."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    from lrgsglib.utils.lrg.infocomm import (
        compute_entropy_observables_from_eigenvalues,
    )

    _, eigvals = get_graph_lspectrum(small_mcg.gr["G"], library="numpy")
    S, C, var, tau = compute_entropy_observables_from_eigenvalues(
        eigenvalues=eigvals, num_nodes=small_mcg.N,
    )

    max_S = np.log(small_mcg.N)
    assert np.all(S >= -0.01), f"Entropy should be >= 0, min={S.min()}"
    assert np.all(S <= max_S + 0.01), f"Entropy should be <= log(N)={max_S:.3f}"


@pytest.mark.physical
def test_specific_heat_non_negative(small_mcg):
    """Specific heat C(tau) should be non-negative (allow small numerical error)."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    from lrgsglib.utils.lrg.infocomm import (
        compute_entropy_observables_from_eigenvalues,
    )

    _, eigvals = get_graph_lspectrum(small_mcg.gr["G"], library="numpy")
    S, C, var, tau = compute_entropy_observables_from_eigenvalues(
        eigenvalues=eigvals, num_nodes=small_mcg.N,
    )
    assert np.all(C >= -0.1), f"Specific heat should be >= 0, min={C.min()}"


@pytest.mark.physical
@pytest.mark.visual
def test_entropy_specific_heat_curves(
    small_mcg, small_lattice_sqr, small_er, output_dir, save_figure,
):
    """Visual: S(tau) and C(tau) for different graph types."""
    import matplotlib.pyplot as plt
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    from lrgsglib.utils.lrg.infocomm import (
        compute_entropy_observables_from_eigenvalues,
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    graphs = {
        "MCG": small_mcg,
        "Lattice2D (sqr)": small_lattice_sqr,
        "Erdos-Renyi": small_er,
    }

    for label, sg in graphs.items():
        _, eigvals = get_graph_lspectrum(sg.gr["G"], library="numpy")
        S, C, var, tau = compute_entropy_observables_from_eigenvalues(
            eigenvalues=eigvals, num_nodes=sg.N,
        )
        ax1.semilogx(tau, S, label=label)
        ax2.semilogx(tau, C, label=label)

    ax1.set_xlabel(r"$\tau$")
    ax1.set_ylabel(r"$S(\tau)$")
    ax1.set_title("Von Neumann Entropy")
    ax1.legend()

    ax2.set_xlabel(r"$\tau$")
    ax2.set_ylabel(r"$C(\tau)$")
    ax2.set_title("Specific Heat")
    ax2.legend()

    fig.tight_layout()
    save_figure(fig, output_dir / "entropy_specific_heat_comparison.png")


# ===================================================================
# Density of states (DOS)
# ===================================================================


@pytest.mark.physical
@pytest.mark.visual
def test_dos_comparison(
    small_mcg, small_lattice_sqr, small_er, output_dir, save_figure,
):
    """Visual: density of states for different graph types."""
    import matplotlib.pyplot as plt
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

    fig, ax = plt.subplots(figsize=(8, 5))

    graphs = {
        "MCG": small_mcg,
        "Lattice2D (sqr)": small_lattice_sqr,
        "Erdos-Renyi": small_er,
    }

    for label, sg in graphs.items():
        _, eigvals = get_graph_lspectrum(sg.gr["G"], library="numpy")
        ax.hist(eigvals, bins=30, alpha=0.5, density=True, label=label)

    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel("Density")
    ax.set_title("Density of States (Laplacian)")
    ax.legend()

    fig.tight_layout()
    save_figure(fig, output_dir / "dos_comparison.png")


# ===================================================================
# Inverse Participation Ratio (IPR)
# ===================================================================


@pytest.mark.physical
def test_ipr_delocalized(small_lattice_sqr):
    """On a regular lattice, eigenvectors should be mostly delocalized (IPR ~ 1/N)."""
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

    L, eigvals = get_graph_lspectrum(
        small_lattice_sqr.gr["G"], library="numpy",
    )
    # Compute eigenvectors
    _, eigvecs = np.linalg.eigh(L)

    N = small_lattice_sqr.N
    # IPR = sum(|psi_i|^4); for delocalized state IPR ~ 1/N
    iprs = np.array([np.sum(np.abs(eigvecs[:, k]) ** 4) for k in range(N)])

    median_ipr = np.median(iprs)
    # Delocalized: IPR should be O(1/N), allow factor of 10
    assert median_ipr < 10.0 / N, (
        f"Median IPR={median_ipr:.4f} too large for delocalized states "
        f"(expected ~ 1/N = {1.0 / N:.4f})"
    )
