"""pb-vs-py agreement on SIGNED graphs for the vector quad (audit gap #1).

The pb kernels receive the signed weight array, but until now every
pb-vs-py comparison ran at pflip=0 — nothing verified the native kernels
on a frustrated substrate. For each model: same signed lattice
(pflip=0.15), several seeds per backend, compare the tail-averaged
equilibrium energy per spin (D-B6: statistical agreement, never
byte-equality — the RNG streams differ by construction).
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D

SIDE = 8
PFLIP = 0.15
STEPS = 200
TAIL = 50
SEEDS = (1, 2, 3, 4)


def make_lattice(tmp_path, seed=7):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=PFLIP,
        engine="nx",
        path_data=tmp_path,
        seed=seed,
    )


def _tail_energy(cls, tmp_path, runlang, T, **kw):
    vals = []
    for sd in SEEDS:
        sg = make_lattice(tmp_path)  # same graph seed = same disorder
        m = cls(
            sg,
            T=T,
            steps=STEPS,
            seed=sd,
            savedisk=False,
            runlang=runlang,
            **kw,
        )
        m.run()
        vals.append(float(np.mean(np.asarray(m.ene)[-TAIL:])))
    return float(np.mean(vals)), float(np.std(vals))


def _assert_backends_agree(cls, tmp_path, T, tol, **kw):
    e_py, s_py = _tail_energy(cls, tmp_path, "py", T, **kw)
    e_pb, s_pb = _tail_energy(cls, tmp_path, "pb", T, **kw)
    scatter = max(s_py, s_pb, 1e-3)
    assert abs(e_py - e_pb) < max(tol, 4.0 * scatter), (
        f"{cls.__name__} signed pb/py disagree: py={e_py:.4f}±{s_py:.4f} "
        f"pb={e_pb:.4f}±{s_pb:.4f}"
    )


def _needs(modname):
    def _check():
        try:
            __import__(modname)
            return True
        except Exception:
            return False

    return pytest.mark.skipif(not _check(), reason=f"{modname} not built")


@_needs("lrgsglib.statsys.PottsModel.ccore._potts_native")
def test_potts_pb_signed_agreement(tmp_path):
    from lrgsglib.statsys.PottsModel import PottsMetropolis

    _assert_backends_agree(PottsMetropolis, tmp_path, T=1.0, tol=0.08, q=3)


@_needs("lrgsglib.statsys.XYModel.ccore._xy_native")
def test_xy_pb_signed_agreement(tmp_path):
    from lrgsglib.statsys.XYModel import XYMetropolis

    _assert_backends_agree(XYMetropolis, tmp_path, T=0.8, tol=0.08)


@_needs("lrgsglib.statsys.HeisenbergModel.ccore._heisenberg_native")
def test_heisenberg_pb_signed_agreement(tmp_path):
    from lrgsglib.statsys.HeisenbergModel import HeisenbergMetropolis

    _assert_backends_agree(HeisenbergMetropolis, tmp_path, T=0.8, tol=0.08)


@_needs("lrgsglib.statsys.MultiSpeciesModel.ccore._multispec_native")
def test_multispec_pb_signed_agreement(tmp_path):
    """pb is identity-M/uniform-q/energy-only, so deselect magn."""
    from lrgsglib.statsys.MultiSpeciesModel import MultiSpeciesMetropolis

    _assert_backends_agree(
        MultiSpeciesMetropolis,
        tmp_path,
        T=1.0,
        tol=0.08,
        species=2,
        q_per_species=3,
        observables=("energy",),
    )
