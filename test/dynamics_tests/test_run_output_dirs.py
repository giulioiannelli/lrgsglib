"""Phase C output contract: one directory per run + cfg.json sidecar.

Pins the interview-approved layout::

    <dynpath>/p=0.1_T=2.27_ns=5_lang=py_s=91731/
        cfg.json  ene.bin  m.bin  [sout.bin]

- dirname tokens: numerics first, key=value, defaults elided,
  lang= and s= always written, out_suffix appended verbatim;
- cfg.json: graph / model / run blocks — the reproducibility ground
  truth (graph seed pins the disorder realization);
- construction stays lazy (no directories before run()).
"""

from __future__ import annotations

import json

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys import IsingModel
from lrgsglib.statsys.PottsModel import PottsMetropolis

SIDE = 8
N = SIDE * SIDE


def make_lattice(tmp_path, pflip=0.1):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=7,
    )


def test_run_dir_name_and_contents(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingModel(sg, T=2.27, steps=5, seed=91731)
    assert [p for p in tmp_path.rglob("*") if p.is_dir()] == []  # lazy
    model.run(tqdm_on=False)

    rundir = model._run_output_dir()
    assert rundir.name == "p=0.1_T=2.27_ns=5_lang=py_s=91731"
    assert rundir.parent == model.dynpath
    assert sorted(p.name for p in rundir.iterdir()) == [
        "cfg.json",
        "ene.bin",
        "m.bin",
    ]
    ene = np.fromfile(rundir / "ene.bin", dtype=np.float64)
    assert ene.shape == (5,)


def test_cfg_sidecar_blocks(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingModel(sg, T=1.5, steps=3, seed=11)
    model.run(tqdm_on=False)
    cfg = json.loads((model._run_output_dir() / "cfg.json").read_text())

    assert set(cfg) == {"graph", "model", "run"}
    # graph block: p names the ensemble, the seed pins the realization
    assert cfg["graph"]["pflip"] == 0.1
    assert cfg["graph"]["seed"] == 7
    assert cfg["graph"]["N"] == N
    assert cfg["graph"]["class"] == "Lattice2DNX"
    assert cfg["model"]["class"] == "IsingMetropolis"
    assert cfg["run"]["seed"] == 11
    assert cfg["run"]["steps"] == 3
    assert cfg["run"]["backend"] == "py"
    assert "date" in cfg["run"] and "lrgsglib" in cfg["run"]


def test_out_suffix_appended_verbatim(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingModel(sg, T=1.5, steps=3, seed=5, out_suffix="run17")
    model.run(tqdm_on=False)
    assert model._run_output_dir().name.endswith("_s=5_run17")


def test_repeat_runs_with_different_seeds_do_not_clobber(tmp_path):
    sg = make_lattice(tmp_path)
    dirs = []
    for seed in (1, 2, 3):
        m = PottsMetropolis(sg, q=3, T=1.0, steps=3, seed=seed)
        m.run(tqdm_on=False)
        dirs.append(m._run_output_dir())
    assert len({d.name for d in dirs}) == 3  # seed token disambiguates
    for d in dirs:
        assert (d / "ene.bin").exists()


def test_sweep_glob_pattern(tmp_path):
    """The pattern from the interview: p=0.1_T=*_ns=3_lang=py_s=*/ene.bin"""
    sg = make_lattice(tmp_path)
    for T in (1.0, 2.0, 3.0):
        IsingModel(sg, T=T, steps=3, seed=9).run(tqdm_on=False)
    hits = sorted(sg.path_ising.glob("p=0.1_T=*_ns=3_lang=py_s=*/ene.bin"))
    assert len(hits) == 3
