"""Distribution-level correctness of cluster moves on a SIGNED graph.

Audit gap #3: cluster tests asserted "runs + energy bookkeeping", never
that Wolff/Swendsen-Wang actually sample the Boltzmann distribution.
Here the substrate is small enough (3x3 periodic lattice, N=9, 512
states) to enumerate exp(-E/T) EXACTLY: the sampled time averages of the
energy and |magnetization| must match the exact ensemble averages within
a few autocorrelation-corrected standard errors — on a frustrated
(pflip>0) graph, where the embedded-Ising bookkeeping is nontrivial.
"""

from __future__ import annotations

import itertools

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.IsingDynamics import IsingMetropolis

SIDE = 3
N = SIDE * SIDE
T = 2.5
STEPS = 6000
BURN = 1000
# generous integrated-autocorrelation allowance for the error bars
TAU = 10.0


def _exact_ensemble(model):
    """Enumerate all 2^N states: exact <e>, <|m|>, and their variances."""
    states = np.array(list(itertools.product((-1, 1), repeat=N)), dtype=np.int8)
    E = np.array([model.total_energy(s) for s in states])  # EXTENSIVE
    logw = -E / T
    w = np.exp(logw - logw.max())
    w /= w.sum()
    e = E / float(N)  # the recorded series is intensive (plan section 3b)
    absm = np.abs(states.mean(axis=1))
    stats = {}
    for name, x in (("e", e), ("absm", absm)):
        mean = float(np.sum(w * x))
        var = float(np.sum(w * x**2) - mean**2)
        stats[name] = (mean, var)
    return stats


@pytest.mark.physical
@pytest.mark.parametrize("move", ["wolff", "sw"])
def test_cluster_moves_sample_boltzmann_on_signed_graph(tmp_path, move):
    sg = Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=0.2,
        engine="nx",
        path_data=tmp_path,
        seed=3,
    )
    model = IsingMetropolis(
        sg, T=T, move=move, steps=STEPS, seed=11, savedisk=False
    )
    exact = _exact_ensemble(model)
    assert np.any(model._nbr_J < 0.0)  # the substrate is truly frustrated

    model.run(tqdm_on=False)
    e_t = np.asarray(model.ene, dtype=float)[BURN:]
    m_t = np.abs(np.asarray(model.magn, dtype=float))[BURN:]
    n = e_t.size

    for name, sampled in (("e", e_t), ("absm", m_t)):
        mean_exact, var_exact = exact[name]
        stderr = np.sqrt(var_exact * TAU / n)
        assert abs(sampled.mean() - mean_exact) < 5.0 * stderr + 1e-12, (
            f"{move}/{name}: sampled {sampled.mean():.5f} vs exact "
            f"{mean_exact:.5f} (allowed ±{5.0 * stderr:.5f})"
        )


@pytest.mark.physical
def test_single_site_reference_samples_boltzmann_too(tmp_path):
    """Anchor: the single-site sampler passes the identical check, so a
    cluster failure above cannot be blamed on the test's tolerances."""
    sg = Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=0.2,
        engine="nx",
        path_data=tmp_path,
        seed=3,
    )
    model = IsingMetropolis(sg, T=T, steps=STEPS, seed=13, savedisk=False)
    exact = _exact_ensemble(model)
    model.run(tqdm_on=False)
    e_t = np.asarray(model.ene, dtype=float)[BURN:]
    mean_exact, var_exact = exact["e"]
    stderr = np.sqrt(var_exact * TAU / e_t.size)
    assert abs(e_t.mean() - mean_exact) < 5.0 * stderr
