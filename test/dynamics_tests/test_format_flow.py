"""Phase-C output contract for the flow trio (Kuramoto / RD / CoupledODE).

Pins the per-run layout on the graph-anchored dynamics tree::

    data/<graph>/<model>/N=.../<run dirname>/{cfg.json, sfin.bin, ...}

- one hand-computed dirname per model (locked token schema: disorder p
  first, physics numerics always written, dt, ns, categoricals elided at
  their defaults, lang= and s= always last);
- cfg.json carries the graph / model / run blocks;
- sfin.bin is always written under ``savedisk`` (N float64);
- parse_run_dirname round-trips the dirname;
- pb backend: ``lang=pb`` in the dirname + py-vs-pb final-state parity
  (deterministic integrators);
- Kuramoto C backend: the C exchange is redirected into the rundir and
  the order-parameter read-back works.

All writes go under ``tmp_path`` (never the repo's ``data/``).
"""

from __future__ import annotations

import json

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys._naming import parse_run_dirname
from lrgsglib.statsys._solver import SolverBackend
from lrgsglib.statsys._solver_engine import is_backend_available
from lrgsglib.statsys.CoupledODEModel import CoupledODEModel
from lrgsglib.statsys.CoupledODEModel.defaults import CODE_SOLVER_NAME
from lrgsglib.statsys.KuramotoModel import KuramotoModel
from lrgsglib.statsys.KuramotoModel.defaults import KURAMOTO_SOLVER_NAME
from lrgsglib.statsys.ReactionDiffusionModel import ReactionDiffusionModel
from lrgsglib.statsys.ReactionDiffusionModel.defaults import RD_SOLVER_NAME

SIDE = 4
N = SIDE * SIDE

_needs_kuramoto_pb = pytest.mark.skipif(
    not is_backend_available(KURAMOTO_SOLVER_NAME, SolverBackend.PB),
    reason="_kuramoto_native not built",
)
_needs_rd_pb = pytest.mark.skipif(
    not is_backend_available(RD_SOLVER_NAME, SolverBackend.PB),
    reason="_rd_native not built",
)
_needs_code_pb = pytest.mark.skipif(
    not is_backend_available(CODE_SOLVER_NAME, SolverBackend.PB),
    reason="_code_native not built",
)
_needs_kuramoto_c = pytest.mark.skipif(
    not is_backend_available(KURAMOTO_SOLVER_NAME, SolverBackend.C),
    reason="KuramotoSimulator binary not built",
)


def make_lattice(tmp_path, pflip=0.1):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=pflip,
        engine="nx",
        path_data=tmp_path,
        seed=7,
    )


def _kuramoto(tmp_path, *, runlang, steps=400, seed=5, pflip=0.1):
    sg = make_lattice(tmp_path, pflip=pflip)
    return KuramotoModel(
        sg,
        coupling=4.0,
        omega=0.0,
        dt=0.05,
        steps=steps,
        save_order_param=True,
        runlang=runlang,
        seed=seed,
    )


def _rd(tmp_path, *, runlang, steps=1000, seed=3, pflip=0.1):
    sg = make_lattice(tmp_path, pflip=pflip)
    return ReactionDiffusionModel(
        sg,
        D=0.1,
        reaction="fisher_kpp",
        reaction_params={"r": 1.0},
        dt=0.01,
        steps=steps,
        runlang=runlang,
        seed=seed,
    )


def _ode(tmp_path, *, runlang, steps=200, seed=7, pflip=0.1):
    sg = make_lattice(tmp_path, pflip=pflip)
    return CoupledODEModel(
        sg,
        coupling_type="product",
        local_type="fitzhugh",
        coupling_strength=0.5,
        dt=0.01,
        steps=steps,
        runlang=runlang,
        seed=seed,
    )


def _check_run_outputs(model, expected_dirname):
    rundir = model._run_output_dir()
    assert rundir.name == expected_dirname
    assert rundir.parent == model.dynpath
    cfg = json.loads((rundir / "cfg.json").read_text())
    assert set(cfg) == {"graph", "model", "run"}
    assert cfg["model"]["class"] == type(model).__name__
    sfin = np.fromfile(rundir / "sfin.bin", dtype=np.float64)
    assert sfin.shape == (N,)
    assert np.allclose(sfin, np.asarray(model.s, dtype=np.float64))
    # dirname round-trip at the string level
    parsed = parse_run_dirname(rundir.name)
    assert parsed["p"] == "0.1"
    assert parsed["lang"] == "py"
    return rundir, parsed


# ======================================================================
# Pinned dirnames (py backend, one hand-computed case per model)
# ======================================================================


def test_kuramoto_pinned_dirname_and_outputs(tmp_path):
    km = _kuramoto(tmp_path, runlang="py")
    km.run(tqdm_on=False)
    rundir, parsed = _check_run_outputs(
        km, "p=0.1_K=4_om=0_dt=0.05_ns=400_lang=py_s=5"
    )
    assert parsed["s"] == "5" and parsed["ns"] == "400"
    # save_order_param=True -> the order-parameter series is persisted
    r = np.fromfile(rundir / "r.bin", dtype=np.float64)
    assert r.shape == (400,)
    assert np.all((r >= 0.0) & (r <= 1.0 + 1e-9))


def test_rd_pinned_dirname_and_outputs(tmp_path):
    rd = _rd(tmp_path, runlang="py")
    rd.run()
    # reaction fisher_kpp (default) and intg rk4 (default) are elided;
    # the reaction params surface as sorted key=value pairs after D.
    _check_run_outputs(rd, "p=0.1_D=0.1_r=1_dt=0.01_ns=1000_lang=py_s=3")


def test_ode_pinned_dirname_and_outputs(tmp_path):
    ode = _ode(tmp_path, runlang="py")
    ode.run()
    # cpl=prod and loc=fitz are non-default -> written (VOCAB-shortened);
    # intg rk4 (default) elided.
    _check_run_outputs(
        ode, "p=0.1_K=0.5_dt=0.01_ns=200_cpl=prod_loc=fitz_lang=py_s=7"
    )


# ======================================================================
# pb backend: lang=pb in the dirname + py-vs-pb final-state parity
# ======================================================================


@_needs_kuramoto_pb
def test_kuramoto_pb_dirname_and_parity(tmp_path):
    # Construct-and-run sequentially: the model seeds the GLOBAL numpy
    # RNG at construction, so interleaving would desync the two ICs.
    mpy = _kuramoto(tmp_path / "py", runlang="py", steps=100)
    mpy.run(tqdm_on=False)
    mpb = _kuramoto(tmp_path / "pb", runlang="pb", steps=100)
    mpb.run(tqdm_on=False)
    assert np.allclose(mpy.sini, mpb.sini)  # same initial condition
    assert parse_run_dirname(mpb._run_output_dir().name)["lang"] == "pb"
    spy = np.fromfile(mpy._run_output_dir() / "sfin.bin", dtype=np.float64)
    spb = np.fromfile(mpb._run_output_dir() / "sfin.bin", dtype=np.float64)
    # wrap-aware phase comparison (angles live in [0, 2*pi))
    dtheta = np.abs((spy - spb + np.pi) % (2 * np.pi) - np.pi)
    assert np.allclose(dtheta, 0.0, atol=1e-8)


@_needs_rd_pb
def test_rd_pb_dirname_and_parity(tmp_path):
    # Unsigned lattice: the native RD kernel's signed-Laplacian convention
    # diverges from the Python one on flipped edges (pre-existing), so the
    # numerical py<->pb bar applies at pflip=0.
    mpy = _rd(tmp_path / "py", runlang="py", steps=100, pflip=0.0)
    mpy.run()
    mpb = _rd(tmp_path / "pb", runlang="pb", steps=100, pflip=0.0)
    mpb.run()
    assert np.allclose(mpy.sini, mpb.sini)  # same initial condition
    assert parse_run_dirname(mpb._run_output_dir().name)["lang"] == "pb"
    spy = np.fromfile(mpy._run_output_dir() / "sfin.bin", dtype=np.float64)
    spb = np.fromfile(mpb._run_output_dir() / "sfin.bin", dtype=np.float64)
    assert np.allclose(spy, spb, atol=1e-9)


@_needs_code_pb
def test_ode_pb_dirname_and_parity(tmp_path):
    # linear coupling matches py<->pb numerically on an UNSIGNED lattice
    # (signed 'linear' degree conventions differ, see test_code_backend).
    def build(sub, runlang):
        sg = make_lattice(tmp_path / sub, pflip=0.0)
        return CoupledODEModel(
            sg,
            coupling_type="linear",
            local_type="fitzhugh",
            coupling_strength=0.5,
            dt=0.01,
            steps=100,
            runlang=runlang,
            seed=7,
        )

    mpy = build("py", "py")
    mpy.run()
    mpb = build("pb", "pb")
    mpb.run()
    assert np.allclose(mpy.sini, mpb.sini)  # same initial condition
    assert parse_run_dirname(mpb._run_output_dir().name)["lang"] == "pb"
    spy = np.fromfile(mpy._run_output_dir() / "sfin.bin", dtype=np.float64)
    spb = np.fromfile(mpb._run_output_dir() / "sfin.bin", dtype=np.float64)
    assert np.allclose(spy, spb, atol=1e-9)


# ======================================================================
# C backend (Kuramoto): rundir-redirected exchange + read-back
# ======================================================================


@_needs_kuramoto_c
def test_kuramoto_c_rundir_exchange_and_readback(tmp_path):
    mc = _kuramoto(tmp_path, runlang="C0", steps=50)
    mc.run()
    rundir = mc._run_output_dir()
    assert parse_run_dirname(rundir.name)["lang"] == "c"
    names = sorted(p.name for p in rundir.iterdir())
    # cfg sidecar + final state + the C-written order-parameter file;
    # the transient exchange files (s_/h_/edgelist) were cleaned up.
    assert "cfg.json" in names
    assert "sfin.bin" in names
    c_order = mc.sg.get_p_fname("r", mc.out_suffix)
    assert str(c_order) in names
    assert not any(n.startswith(("s_", "h_")) for n in names)
    # read-back populated the in-RAM series (and r.bin was persisted)
    assert len(mc.order_params) == 50
    r = np.fromfile(rundir / "r.bin", dtype=np.float64)
    assert r.shape == (50,)
    assert np.allclose(r, np.asarray(mc.order_params))
    sfin = np.fromfile(rundir / "sfin.bin", dtype=np.float64)
    assert sfin.shape == (N,) and np.all(np.isfinite(sfin))
