"""Signed-graph behavior of the flow models (audit gap: zero signed tests).

Kuramoto couples through the SIGNED adjacency (a negative edge is
anti-phase coupling); ReactionDiffusion diffuses on the SIGNED
Laplacian. Neither had a single pflip>0 test — these pin that the edge
signs actually reach the dynamics, with an unsigned control beside each
signed case.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D

SIDE = 6
STEPS = 400
DT = 0.05


def make_lattice(tmp_path, pflip):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=11,
    )


# ======================================================================
# Kuramoto: signed adjacency
# ======================================================================


def _kuramoto_final_r(tmp_path, pflip, seed=5):
    from lrgsglib.statsys.KuramotoModel import KuramotoModel

    sg = make_lattice(tmp_path, pflip)
    km = KuramotoModel(
        sg,
        coupling=4.0,
        omega=0.0,
        dt=DT,
        steps=STEPS,
        runlang="py",
        seed=seed,
    )
    km.run(tqdm_on=False)
    return km.order_parameter()


@pytest.mark.physical
def test_kuramoto_ferromagnetic_synchronizes(tmp_path):
    """Unsigned control: strong coupling, identical frequencies -> r ~ 1."""
    assert _kuramoto_final_r(tmp_path / "pos", pflip=0.0) > 0.9


@pytest.mark.physical
def test_kuramoto_all_negative_edges_antisynchronize(tmp_path):
    """pflip=1 on the bipartite square lattice: every edge is anti-phase
    coupling, so the two sublattices lock pi apart and the GLOBAL order
    parameter collapses. If the kernel ignored edge signs this would
    synchronize exactly like the control above."""
    assert _kuramoto_final_r(tmp_path / "neg", pflip=1.0) < 0.4


@pytest.mark.physical
def test_kuramoto_negative_edge_pair_locks_antiphase(tmp_path):
    """Smallest possible case: check the sign enters the coupling term."""
    from lrgsglib.statsys.KuramotoModel import KuramotoModel

    sg = make_lattice(tmp_path, pflip=1.0)
    km = KuramotoModel(
        sg,
        coupling=4.0,
        omega=0.0,
        dt=DT,
        steps=STEPS,
        runlang="py",
        seed=5,
    )
    km.run(tqdm_on=False)
    theta = np.asarray(km.s, dtype=float)
    # nearest-neighbour phase differences concentrate near pi (mod 2pi)
    u = theta.reshape(SIDE, SIDE)
    diffs = np.concatenate(
        [
            (u - np.roll(u, 1, axis=0)).ravel(),
            (u - np.roll(u, 1, axis=1)).ravel(),
        ]
    )
    anti = np.abs(np.cos(diffs))  # |cos| ~ 1 when locked at 0 or pi
    assert np.mean(np.cos(diffs)) < -0.7  # locked near pi, not near 0
    assert np.mean(anti) > 0.8


# ======================================================================
# ReactionDiffusion: signed Laplacian
# ======================================================================


def _rd_run(tmp_path, pflip, seed=3):
    from lrgsglib.statsys.ReactionDiffusionModel import ReactionDiffusionModel

    sg = make_lattice(tmp_path, pflip)
    rd = ReactionDiffusionModel(
        sg,
        D=0.1,
        reaction="none",
        dt=0.01,
        steps=STEPS,
        runlang="py",
        seed=seed,
        ic="random",
    )
    rd.init_rd_dynamics()  # realize the ic now so u0 is observable
    u0 = np.asarray(rd.s, dtype=float).copy()
    rd.run(tqdm_on=False)
    return u0, np.asarray(rd.s, dtype=float)


@pytest.mark.physical
def test_rd_unsigned_diffusion_conserves_mass_and_smooths(tmp_path):
    """Unsigned control, no reaction: heat flow conserves the total mass
    (uniform mode = Laplacian kernel) and contracts spatial variance."""
    u0, u1 = _rd_run(tmp_path / "pos", pflip=0.0)
    assert u1.mean() == pytest.approx(u0.mean(), rel=0.05)
    assert np.std(u1) < 0.6 * np.std(u0)


@pytest.mark.physical
def test_rd_signed_laplacian_decays_uniform_mass(tmp_path):
    """pflip=1 on the bipartite lattice: the signed Laplacian is the
    signless Laplacian, whose spectrum is strictly positive on the
    uniform mode — so the total mass now DECAYS, the exact opposite of
    the unsigned control above. (The staggered kernel mode is not
    exactly conserved here because the stepper clamps concentrations to
    [0, inf) each step — a physical nonlinearity.) A sign-blind kernel
    would conserve the mass like the control."""
    u0, u1 = _rd_run(tmp_path / "neg", pflip=1.0)
    assert abs(u1.mean()) < 0.2 * abs(u0.mean())  # uniform mass decays
    assert np.linalg.norm(u1) < 0.3 * np.linalg.norm(u0)  # global damping
