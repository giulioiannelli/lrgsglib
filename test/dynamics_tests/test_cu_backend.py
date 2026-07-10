"""cu (CuPy/GPU) backend physics tests for the spin quad.

The cu backend is the SAME VectorSyncEngine as np, compiled against the
CuPy array module with the substrate staged on device (``_vec_bind``)
and mirrored back each sweep for the host-side observables. These tests
require a working GPU (skipped otherwise) and pin:

- py-vs-cu statistical agreement of the equilibrium energy on a SIGNED
  lattice, per model (mirrors ``test_np_backend``);
- the incremental-energy invariant after a device run;
- the async hard-refusal (sync is the only vectorizable schedule).
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")
cupy = pytest.importorskip("cupy")

try:  # a driverless machine imports cupy fine but cannot allocate
    cupy.zeros(1) + 1
except Exception:  # pragma: no cover - depends on local GPU
    pytest.skip("no usable GPU device", allow_module_level=True)

from lrgsglib.graphs import Lattice2D

SIDE = 8
PFLIP = 0.1
STEPS = 200
TAIL = 80
SEEDS = (11, 12, 13, 14)
TOL = 0.08


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
    if name == "ising":
        from lrgsglib.statsys.IsingDynamics import IsingMetropolis

        return IsingMetropolis(sg=sg, T=1.5, **kw)
    if name == "potts":
        from lrgsglib.statsys.PottsModel import PottsMetropolis

        return PottsMetropolis(sg=sg, q=3, T=1.0, **kw)
    if name == "xy":
        from lrgsglib.statsys.XYModel import XYMetropolis

        return XYMetropolis(sg=sg, T=0.8, **kw)
    from lrgsglib.statsys.HeisenbergModel import HeisenbergMetropolis

    return HeisenbergMetropolis(sg=sg, T=0.8, **kw)


MODELS = ("heisenberg", "ising", "potts", "xy")


def _tail_energies(name, tmp_path, runlang):
    tails = []
    for seed in SEEDS:
        sg = make_lattice(tmp_path / f"{runlang}{seed}")
        m = _make_model(
            name,
            sg,
            upd_mode="sync",
            steps=STEPS,
            seed=seed,
            savedisk=False,
            runlang=runlang,
        )
        m.run(tqdm_on=False)
        assert abs(m._e_running - m.total_energy()) < 1e-6
        tails.append(float(np.mean(np.asarray(m.ene)[-TAIL:])))
    return np.array(tails)


@pytest.mark.physical
@pytest.mark.parametrize("name", sorted(MODELS))
def test_cu_matches_py_sync_energy(name, tmp_path):
    """py-sync and cu-sync sample the same Gibbs measure on a signed
    lattice: tail energies agree within the seed-to-seed scatter."""
    e_py = _tail_energies(name, tmp_path, "py")
    e_cu = _tail_energies(name, tmp_path, "cu")
    scatter = max(e_py.std(), e_cu.std(), 1e-3)
    diff = abs(e_py.mean() - e_cu.mean())
    assert diff < max(TOL, 4.0 * scatter), (
        f"{name} cu/py sync disagree: py={e_py.mean():.4f}±{e_py.std():.4f} "
        f"cu={e_cu.mean():.4f}±{e_cu.std():.4f}"
    )


@pytest.mark.parametrize("name", sorted(MODELS))
def test_cu_default_async_raises(name, tmp_path):
    """runlang='cu' with the default (async) schedule is a hard
    capability error — the vectorized backend implements sync only."""
    sg = make_lattice(tmp_path)
    m = _make_model(name, sg, steps=STEPS, seed=3, savedisk=False, runlang="cu")
    with pytest.raises(NotImplementedError):
        m.run(tqdm_on=False)


def test_cu_lang_token_in_run_dirname(tmp_path):
    """A persisted cu run names its per-run directory with lang=cu."""
    from lrgsglib.statsys.IsingDynamics import IsingMetropolis

    sg = make_lattice(tmp_path)
    m = IsingMetropolis(
        sg=sg, T=1.5, upd_mode="sync", steps=5, seed=3, runlang="cu"
    )
    m.run(tqdm_on=False)
    rundirs = [d for d in sg.path_ising.iterdir() if d.is_dir()]
    assert len(rundirs) == 1
    assert "lang=cu" in rundirs[0].name
