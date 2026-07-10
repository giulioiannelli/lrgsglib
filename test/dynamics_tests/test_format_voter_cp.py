"""Phase-C output contract for VoterModel and ContactProcess.

Pins the per-run directory layout for the two file-exchanging models::

    <dynpath>/p=0.1_ns=100_lang=py_s=42/{m.bin, sout.bin, cldist.npz, cfg.json}

- dirname tokens: disorder p first, physics numerics, run length ns,
  categorical axes (defaults elided), ``lang=`` / ``s=`` always last;
- in-process backends write the SHORT canonical names, the C subprocess
  keeps its own sprintf-composed names INSIDE the same rundir (redirected
  via the composed syshape argument);
- cfg.json sidecar: graph / model / run blocks for every backend;
- ``parse_run_dirname`` round-trips every produced dirname.

All output goes to pytest ``tmp_path`` — NEVER the real ``data/`` tree.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys._naming import parse_run_dirname
from lrgsglib.statsys.ContactProcess import ContactProcessEI, ContactProcessSIR
from lrgsglib.statsys.VoterModel import VoterModel

SIDE = 8
N = SIDE * SIDE
GRAPH_SEED = 7
PFLIP = 0.1


def make_lattice(tmp_path):
    return Lattice2D(
        side1=SIDE,
        geo="sqr",
        pflip=PFLIP,
        engine="nx",
        path_data=tmp_path,
        seed=GRAPH_SEED,
    )


def _voter_c_binary_missing() -> bool:
    return not (VoterModel._c_bin_dir / "VoterSimulator").exists()


def _cp_ei_binary_missing() -> bool:
    return not (ContactProcessEI._c_bin_dir / "ContactProcessEI").exists()


# ===================================================================
# VoterModel
# ===================================================================


def test_voter_py_run_dir_name_and_contents(tmp_path):
    sg = make_lattice(tmp_path)
    v = VoterModel(sg, steps=100, runlang="py", seed=42, savemagn=True)
    assert [p for p in tmp_path.rglob("*") if p.is_dir()] == []  # lazy
    v.run(tqdm_on=False)

    rundir = v._run_output_dir()
    assert rundir.name == "p=0.1_ns=100_lang=py_s=42"
    assert rundir.parent == Path(v.dynpath)
    assert sorted(p.name for p in rundir.iterdir()) == ["cfg.json", "m.bin"]
    magn = np.fromfile(rundir / "m.bin", dtype=np.float64)
    assert magn.shape == (100,)

    cfg = json.loads((rundir / "cfg.json").read_text())
    assert set(cfg) == {"graph", "model", "run"}
    assert cfg["graph"]["pflip"] == PFLIP
    assert cfg["graph"]["seed"] == GRAPH_SEED
    assert cfg["graph"]["N"] == N
    assert cfg["model"]["class"] == "VoterModel"
    assert cfg["model"]["rule"] == "linear"
    assert cfg["model"]["upd_mode"] == "asynchronous"
    assert cfg["model"]["q"] == 2
    assert cfg["run"]["steps"] == 100
    assert cfg["run"]["seed"] == 42
    assert cfg["run"]["backend"] == "py"


def test_voter_qvoter_dirname_pins_rule_tokens(tmp_path):
    sg = make_lattice(tmp_path)
    v = VoterModel(
        sg,
        steps=5,
        runlang="py",
        seed=42,
        savemagn=True,
        rule="qvoter",
        q=3,
        eps=0.05,
    )
    v.run(tqdm_on=False)
    rundir = v._run_output_dir()
    assert rundir.name == "p=0.1_q=3_eps=0.05_ns=5_rule=qv_lang=py_s=42"
    assert (rundir / "m.bin").exists()


def test_voter_nonlinear_dirname_pins_alpha(tmp_path):
    sg = make_lattice(tmp_path)
    v = VoterModel(
        sg, steps=5, runlang="py", seed=42, rule="nonlinear", alpha=1.5
    )
    assert (
        v._run_output_dir().name
        == "p=0.1_eps=0_alpha=1.5_ns=5_rule=nlin_lang=py_s=42"
    )


@pytest.mark.integration
def test_voter_c_run_dir_and_readback(tmp_path):
    if _voter_c_binary_missing():
        pytest.skip("VoterSimulator binary not built")
    sg = make_lattice(tmp_path)
    steps = 100
    v = VoterModel(
        sg,
        steps=steps,
        runlang="C0",
        seed=42,
        savemagn=True,
        absorbing_check=False,
    )
    v.run(tqdm_on=False)

    rundir = v._run_output_dir()
    assert rundir.name == "p=0.1_ns=100_lang=c_s=42"  # lang= C backend
    assert rundir.exists()
    # The C binary keeps its own sprintf-composed filename INSIDE the rundir
    # and the read-back joins the exact same path.
    assert v.magn_path == rundir / "m_p=0.1.bin"
    assert v.magn_path.exists()
    assert len(v.magn) == steps
    assert (rundir / "cfg.json").exists()
    # Transient C inputs are cleaned; the graph-side edgelist mirror too.
    assert not (rundir / "s_p=0.1.bin").exists()
    assert not (sg.path_graph / rundir.name / "edgelist_p=0.1.bin").exists()


# ===================================================================
# ContactProcess SIR (py / pb)
# ===================================================================


def _sir(tmp_path, runlang, seed=42, steps=50):
    sg = make_lattice(tmp_path)
    return ContactProcessSIR(
        sg,
        mu=0.5,
        beta=1.0,
        runlang=runlang,
        seed=seed,
        steps=steps,
        savedensity=True,
        ic="uniform",
    )


def test_cp_sir_py_run_dir_name_and_contents(tmp_path):
    m = _sir(tmp_path, "py")
    m.run(tqdm_on=False)
    rundir = m._run_output_dir()
    assert rundir.name == "p=0.1_mu=0.5_beta=1_ns=50_lang=py_s=42"
    assert sorted(p.name for p in rundir.iterdir()) == ["cfg.json", "rho.bin"]
    rho = np.fromfile(rundir / "rho.bin", dtype=np.float64)
    assert rho.shape == (50,)
    cfg = json.loads((rundir / "cfg.json").read_text())
    assert cfg["model"]["class"] == "ContactProcessSIR"
    assert cfg["model"]["mu"] == 0.5
    assert cfg["model"]["beta"] == 1.0


def test_cp_sir_pb_same_dirname_modulo_lang(tmp_path):
    pytest.importorskip(
        "lrgsglib.statsys.ContactProcess.ccore._cp_native",
        reason="_cp_native not built",
    )
    py = _sir(tmp_path / "py", "py")
    py.run(tqdm_on=False)
    pb = _sir(tmp_path / "pb", "pb")
    pb.run(tqdm_on=False)

    name_py = py._run_output_dir().name
    name_pb = pb._run_output_dir().name
    assert name_pb == "p=0.1_mu=0.5_beta=1_ns=50_lang=pb_s=42"
    assert name_pb.replace("lang=pb", "lang=py") == name_py
    rho = np.fromfile(pb._run_output_dir() / "rho.bin", dtype=np.float64)
    assert rho.shape == (50,)


# ===================================================================
# ContactProcess EI (C subprocess, density mode)
# ===================================================================


@pytest.mark.integration
def test_cp_ei_c_run_dir_and_density_readback(tmp_path):
    if _cp_ei_binary_missing():
        pytest.skip("ContactProcessEI binary not built")
    sg = make_lattice(tmp_path)
    nls = 20
    m = ContactProcessEI(
        sg, gamma=0.8, runlang="C1D", steps=50, seed=3, num_log_samples=nls
    )
    m.run(verbose=False)

    rundir = m._run_output_dir()
    assert rundir.name == "p=0.1_gam=0.8_ns=50_nls=20_lang=c_s=3"
    # The binary's own density filename is preserved inside the rundir …
    assert (rundir / "rho_p=0.1_gamma=0.8.bin").exists()
    assert (rundir / "cfg.json").exists()
    # … and the read-back populates self.density with num_log_samples points.
    rho = np.asarray(m.density, dtype=np.float64)
    assert rho.size == nls
    assert np.all((rho >= 0.0) & (rho <= 1.0))


# ===================================================================
# parse_run_dirname round-trips every schema
# ===================================================================


def test_parse_run_dirname_round_trip(tmp_path):
    sg = make_lattice(tmp_path)
    voter = VoterModel(
        sg, steps=7, runlang="py", seed=11, rule="qvoter", q=4, eps=0.25
    )
    sir = ContactProcessSIR(
        sg, mu=0.5, beta=2.0, runlang="py", seed=11, steps=7
    )
    ei = ContactProcessEI(
        sg,
        gamma=0.8,
        runlang="py",
        seed=11,
        steps=7,
        output_mode="density",
        num_log_samples=30,
        activation="relu",
        update_mode="cached",
    )
    expectations = {
        voter: {
            "p": "0.1",
            "q": "4",
            "eps": "0.25",
            "ns": "7",
            "rule": "qv",
            "lang": "py",
            "s": "11",
        },
        sir: {
            "p": "0.1",
            "mu": "0.5",
            "beta": "2",
            "ns": "7",
            "lang": "py",
            "s": "11",
        },
        ei: {
            "p": "0.1",
            "gam": "0.8",
            "ns": "7",
            "nls": "30",
            "act": "relu",
            "upd": "cached",
            "lang": "py",
            "s": "11",
        },
    }
    for model, want in expectations.items():
        assert parse_run_dirname(model._run_output_dir().name) == want
