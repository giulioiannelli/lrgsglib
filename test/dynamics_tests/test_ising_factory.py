"""``IsingModel(sg, scheme=..., **params)`` — the single scheme-dispatch front door.

The factory must select the scheme class explicitly (never infer it from
which kwargs are present), forward parameters verbatim, accept the short
``sch=``-token vocabulary and the long spellings, and hard-error on an
unknown scheme naming the accepted values.
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys import IsingModel
from lrgsglib.statsys.IsingDynamics import (
    IsingCEM,
    IsingMetropolis,
    IsingParallelTempering,
    IsingSimulatedAnnealing,
)

SIDE = 8
SEED = 42


def make_lattice(tmp_path, pflip=0.1):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=7,
    )


def test_default_scheme_is_metropolis(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingModel(sg, T=1.5, steps=3, seed=SEED, savedisk=False)
    assert type(model) is IsingMetropolis
    assert model.T == 1.5
    assert model.steps == 3


@pytest.mark.parametrize(
    "scheme, cls, params",
    [
        ("metropolis", IsingMetropolis, {"T": 2.0}),
        ("met", IsingMetropolis, {"T": 2.0}),
        ("sa", IsingSimulatedAnnealing, {"T_init": 5.0, "T_final": 0.5}),
        (
            "simulated_annealing",
            IsingSimulatedAnnealing,
            {"T_init": 5.0, "T_final": 0.5},
        ),
        (
            "pt",
            IsingParallelTempering,
            {"n_replicas": 4, "T_min": 1.0, "T_max": 3.0},
        ),
        (
            "parallel_tempering",
            IsingParallelTempering,
            {"n_replicas": 4, "T_min": 1.0, "T_max": 3.0},
        ),
        ("cem", IsingCEM, {"restarts": 2, "n_iter": 3}),
        ("cross_entropy", IsingCEM, {"restarts": 2, "n_iter": 3}),
    ],
)
def test_scheme_dispatch_and_aliases(tmp_path, scheme, cls, params):
    sg = make_lattice(tmp_path)
    model = IsingModel(sg, scheme=scheme, seed=SEED, savedisk=False, **params)
    assert type(model) is cls


def test_scheme_is_case_insensitive(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingModel(
        sg,
        scheme="PT",
        n_replicas=4,
        T_min=1.0,
        T_max=3.0,
        seed=SEED,
        savedisk=False,
    )
    assert type(model) is IsingParallelTempering


def test_params_forwarded_verbatim(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingModel(
        sg,
        scheme="sa",
        T_init=8.0,
        T_final=0.25,
        n_temperatures=6,
        steps_per_T=2,
        seed=SEED,
        savedisk=False,
    )
    assert model.T_schedule[0] == pytest.approx(8.0)
    assert len(model.T_schedule) == 6
    assert np.all(np.diff(model.T_schedule) < 0.0)


def test_unknown_scheme_raises_naming_options(tmp_path):
    sg = make_lattice(tmp_path)
    with pytest.raises(ValueError, match="unknown Ising scheme.*'met'"):
        IsingModel(sg, scheme="wolff", T=1.0)


def test_factory_instance_runs(tmp_path):
    sg = make_lattice(tmp_path)
    model = IsingModel(sg, T=1.5, steps=3, seed=SEED, savedisk=False)
    model.run(tqdm_on=False)
    assert len(model.ene) == 3
