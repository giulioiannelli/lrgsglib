"""Cross-model output-format contract (Phase C closing sweep).

Every new-style statsys model, run through the shared ``run()`` front
door on a small SIGNED lattice, must land its outputs in the canonical
per-run layout::

    data/<graph>/<model>/N=.../<token dirname>/{<obs>.bin, cfg.json}

with a machine-parseable dirname (``key=value`` tokens; ``p``/``ns``/
``lang``/``s`` always present) and a three-block ``cfg.json`` sidecar.
Exact per-model dirnames are pinned in the per-family format tests
(``test_run_output_dirs``, ``test_format_voter_cp``,
``test_format_flow``); this sweep pins that NO model escapes the
contract.
"""

from __future__ import annotations

import json

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys._naming import parse_run_dirname

SIDE = 4
PFLIP = 0.1
STEPS = 4
SEED = 42
GRAPH_SEED = 7


def _lattice(tmp_path, sub):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=PFLIP,
        engine="nx",
        path_data=tmp_path / sub,
        seed=GRAPH_SEED,
    )


def _ising(sg):
    from lrgsglib.statsys import IsingModel

    return IsingModel(sg, scheme="met", T=2.0, steps=STEPS, seed=SEED)


def _potts(sg):
    from lrgsglib.statsys.PottsModel import PottsMetropolis

    return PottsMetropolis(sg=sg, q=3, T=1.0, steps=STEPS, seed=SEED)


def _xy(sg):
    from lrgsglib.statsys.XYModel import XYMetropolis

    return XYMetropolis(sg=sg, T=1.0, steps=STEPS, seed=SEED)


def _heisenberg(sg):
    from lrgsglib.statsys.HeisenbergModel import HeisenbergMetropolis

    return HeisenbergMetropolis(sg=sg, T=1.0, steps=STEPS, seed=SEED)


def _multispecies(sg):
    from lrgsglib.statsys.MultiSpeciesModel import MultiSpeciesMetropolis

    return MultiSpeciesMetropolis(
        sg,
        species=2,
        q_per_species=3,
        T=1.0,
        steps=STEPS,
        seed=SEED,
    )


def _voter(sg):
    from lrgsglib.statsys.VoterModel import VoterModel

    return VoterModel(sg, steps=STEPS, runlang="py", seed=SEED, savemagn=True)


def _cp_sir(sg):
    from lrgsglib.statsys.ContactProcess import ContactProcessSIR

    return ContactProcessSIR(
        sg,
        mu=0.5,
        steps=STEPS,
        runlang="py",
        seed=SEED,
        savedensity=True,
    )


def _kuramoto(sg):
    from lrgsglib.statsys.KuramotoModel import KuramotoModel

    return KuramotoModel(
        sg,
        coupling=1.0,
        dt=0.01,
        steps=STEPS,
        runlang="py",
        seed=SEED,
        save_order_param=True,
    )


def _rd(sg):
    from lrgsglib.statsys.ReactionDiffusionModel import ReactionDiffusionModel

    return ReactionDiffusionModel(
        sg,
        D=0.1,
        reaction="none",
        dt=0.01,
        steps=STEPS,
        runlang="py",
        seed=SEED,
    )


def _coupled_ode(sg):
    from lrgsglib.statsys.CoupledODEModel import CoupledODEModel

    return CoupledODEModel(sg, dt=0.01, steps=STEPS, runlang="py", seed=SEED)


CASES = [
    ("ising", "path_ising", _ising),
    ("potts", "path_potts", _potts),
    ("xy", "path_xy", _xy),
    ("heisenberg", "path_heisenberg", _heisenberg),
    ("multi_species", "path_multi_species", _multispecies),
    ("voter", "path_voter", _voter),
    ("cp_sir", "path_cntct", _cp_sir),
    ("kuramoto", "path_kuramoto", _kuramoto),
    ("reaction_diffusion", "path_reaction_diffusion", _rd),
    ("coupled_ode", "path_coupled_ode", _coupled_ode),
]


@pytest.mark.parametrize(
    "name,path_attr,build", CASES, ids=[c[0] for c in CASES]
)
def test_run_lands_in_canonical_per_run_layout(
    tmp_path, name, path_attr, build
):
    sg = _lattice(tmp_path, name)
    model = build(sg)
    model.run(tqdm_on=False)

    dynpath = getattr(sg, path_attr)
    rundirs = [
        d for d in dynpath.iterdir() if d.is_dir() and (d / "cfg.json").exists()
    ]
    assert len(rundirs) == 1, (
        f"{name}: expected exactly one per-run dir with cfg.json under "
        f"{dynpath}, found {sorted(p.name for p in dynpath.iterdir())}"
    )
    rundir = rundirs[0]

    # Dirname is machine-parseable and carries the four universal tokens.
    parsed = parse_run_dirname(rundir.name)
    assert parsed["p"] == f"{PFLIP:.3g}"
    assert parsed["ns"] == str(STEPS)
    assert parsed["lang"] == "py"
    assert parsed["s"] == str(SEED)

    # cfg.json is the three-block reproducibility sidecar.
    cfg = json.loads((rundir / "cfg.json").read_text())
    assert set(cfg) == {"graph", "model", "run"}
    assert cfg["graph"]["pflip"] == PFLIP
    assert cfg["graph"]["seed"] == GRAPH_SEED
    assert cfg["model"]["class"] == type(model).__name__
    assert cfg["run"]["steps"] == STEPS
    assert cfg["run"]["seed"] == SEED

    # At least one observable artifact accompanies the sidecar.
    payloads = [p for p in rundir.iterdir() if p.suffix != ".json"]
    assert payloads, f"{name}: rundir holds only cfg.json"
    assert all(p.stat().st_size > 0 for p in payloads)
