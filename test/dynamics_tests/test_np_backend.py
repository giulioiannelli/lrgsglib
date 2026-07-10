"""Vectorized np (VectorSyncEngine) backend tests for the spin quad.

For each of Ising / Potts / XY / Heisenberg (new scheme classes):

- np sync tail energy statistically matches py sync on a signed lattice
  (same tolerance pattern as ``test_pb_signed``: several seeds, tail mean,
  4x scatter tolerance — D-B6, statistical agreement, never byte-equality);
- the vectorized ``_vec_delta_E`` hook matches the scalar ``delta_E``
  site-by-site on a fixed state and fixed proposals;
- ``runlang='np'`` with the default (async) update mode is a hard
  capability error (invariant #3: never a silent fallback);
- the incremental-energy invariant holds after an np run.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys._vectorized import greedy_color_classes
from lrgsglib.statsys.HeisenbergModel import HeisenbergMetropolis
from lrgsglib.statsys.IsingDynamics import IsingMetropolis
from lrgsglib.statsys.PottsModel import PottsMetropolis
from lrgsglib.statsys.XYModel import XYMetropolis

SIDE = 8
PFLIP = 0.1
STEPS = 200
TAIL = 50
SEEDS = (1, 2, 3, 4)
TOL = 0.08
VEC_SEED = 12345

#: model key -> (class, physics kwargs)
MODELS = {
    "ising": (IsingMetropolis, {"T": 1.5}),
    "potts": (PottsMetropolis, {"T": 1.0, "q": 3}),
    "xy": (XYMetropolis, {"T": 0.8}),
    "heisenberg": (HeisenbergMetropolis, {"T": 0.8}),
}


def make_lattice(tmp_path, seed=7):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=PFLIP,
        engine="nx",
        path_data=tmp_path,
        seed=seed,
    )


def _make_model(name, sg, **kw):
    cls, physics = MODELS[name]
    return cls(sg, **{**physics, **kw})


def _tail_energy(name, tmp_path, runlang):
    vals = []
    for sd in SEEDS:
        sg = make_lattice(tmp_path)  # same graph seed = same disorder
        m = _make_model(
            name,
            sg,
            upd_mode="sync",
            steps=STEPS,
            seed=sd,
            savedisk=False,
            runlang=runlang,
        )
        m.run()
        vals.append(float(np.mean(np.asarray(m.ene)[-TAIL:])))
    return float(np.mean(vals)), float(np.std(vals))


@pytest.mark.parametrize("name", sorted(MODELS))
def test_np_sync_matches_py_sync(name, tmp_path):
    e_py, s_py = _tail_energy(name, tmp_path, "py")
    e_np, s_np = _tail_energy(name, tmp_path, "np")
    scatter = max(s_py, s_np, 1e-3)
    assert abs(e_py - e_np) < max(TOL, 4.0 * scatter), (
        f"{name} np/py sync disagree: py={e_py:.4f}±{s_py:.4f} "
        f"np={e_np:.4f}±{s_np:.4f}"
    )


@pytest.mark.parametrize("name", sorted(MODELS))
def test_vec_delta_E_matches_scalar(name, tmp_path):
    """The vectorized ΔE hook must reproduce the scalar substrate ΔE
    site-by-site (same state, same proposals — the physics is identical)."""
    sg = make_lattice(tmp_path)
    m = _make_model(
        name,
        sg,
        upd_mode="sync",
        steps=STEPS,
        seed=3,
        savedisk=False,
        runlang="np",
    )
    np.random.seed(VEC_SEED)
    m.init_dynamics()
    m._vec_bind(np)  # stage the substrate arrays (np: zero-copy aliases)
    classes = greedy_color_classes(m._nbr_ptr, m._nbr_idx, m.N)
    cls = classes[0]
    proposals = m._vec_propose(cls, np, np.random)
    vec = np.asarray(m._vec_delta_E(cls, proposals, np), dtype=np.float64)
    scalar = np.array(
        [m.delta_E(int(nd), proposals[k]) for k, nd in enumerate(cls)]
    )
    assert vec.shape == scalar.shape
    assert np.allclose(vec, scalar, rtol=0.0, atol=1e-9), (
        f"{name} _vec_delta_E deviates from scalar delta_E by "
        f"{np.abs(vec - scalar).max():.3e}"
    )


@pytest.mark.parametrize("name", sorted(MODELS))
def test_np_default_async_raises(name, tmp_path):
    """runlang='np' with the default (async) schedule is a hard capability
    error — the vectorized backend implements sync only."""
    sg = make_lattice(tmp_path)
    m = _make_model(name, sg, steps=STEPS, seed=1, savedisk=False, runlang="np")
    with pytest.raises(NotImplementedError):
        m.run()


@pytest.mark.parametrize("name", sorted(MODELS))
def test_np_incremental_energy_invariant(name, tmp_path):
    """After an np run the incrementally tracked energy equals a full
    recompute (ΔE ↔ E consistency, plan §3b)."""
    sg = make_lattice(tmp_path)
    m = _make_model(
        name,
        sg,
        upd_mode="sync",
        steps=STEPS,
        seed=2,
        savedisk=False,
        runlang="np",
    )
    m.run()
    assert abs(m._e_running - m.total_energy()) < 1e-6
