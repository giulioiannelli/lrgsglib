"""
Regression tests for the 2026-09-30 legacy-disorder hotfix.

- B1 double flip: since 844192d the constructors realize the disorder, so a
  kernel that flips ``nwDict[cell]`` afterwards must build with
  ``disorder=None``. Otherwise ``rand`` comes out all-positive and the
  structured cells get a random overlay on top of their pattern.
- B2: the GT kernels must not pin ``seed=42``. Independent calls must draw
  independent realizations.
- B5: concurrent ``write_cfg_sidecar`` writers must not collide.
- B7: run dirnames and sbatch argv keep distinct grid points distinct.
- B9: temp names must not follow the stdlib ``random`` stream, which graph
  constructors reseed.
"""

from __future__ import annotations

import json
import random
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT_DIR = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT_DIR / "src"

# Add src/ to path so `from kernels... import ...` works
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

pytestmark = [pytest.mark.integration]

SIDE = 16
SIDE_3D = 6
PFLIP = 0.1
ER_N, ER_P = 200, 0.05
N_WRITERS = 32


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _key(u, v) -> tuple[int, int]:
    return (int(u), int(v)) if int(u) <= int(v) else (int(v), int(u))


def _negative_edges(sg, on_g: str = "G") -> set[tuple[int, int]]:
    edges = sg.gr[on_g].edges(data="weight")
    return {_key(u, v) for u, v, w in edges if w < 0}


def _assert_pattern_only(sg, cell: str) -> None:
    """The realized negatives are exactly the ``cell`` pattern: not empty
    (``rand`` toggled back) and no random overlay (structured cells)."""
    negatives = _negative_edges(sg)
    assert negatives, f"{cell}: no negative edge left (double flip)"
    assert negatives == {_key(u, v) for u, v in sg.nwDict[cell]["G"]}


def _edge_signature(sg) -> np.ndarray:
    u, v, w = sg.get_edge_arrays()
    return np.stack([u, v, w])


# ---------------------------------------------------------------------------
# B1 - every kernel realizes the disorder exactly once
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("cell", ["randXERR", "randZERR"])
def test_l2d_dict_args_structured_cell_is_pattern_only(cell, tmp_path):
    from kernels.L2D import initialize_l2d_dict_args
    from lrgsglib.graphs.nx import Lattice2DNX

    args = SimpleNamespace(
        L=SIDE, geometry="squared", workdir=None, p=PFLIP, cell_type=cell
    )
    lattice = Lattice2DNX(
        **initialize_l2d_dict_args(args), path_data=tmp_path, path_plot=tmp_path
    )
    lattice.flip_sel_edges(lattice.nwDict[cell]["G"])
    _assert_pattern_only(lattice, cell)


@pytest.mark.parametrize("cell", ["rand", "randXERR"])
def test_l2d_transcluster_lattice_flips_once(cell, tmp_path):
    from kernels.L2D_TransCluster import _make_lattice
    from kernels.TransClusters import build_geometry_selector

    args = SimpleNamespace(L=SIDE, geometry="squared", p=PFLIP)
    lattice = _make_lattice(args, "nx", path_data=tmp_path, path_plot=tmp_path)
    lattice.flip_sel_edges(build_geometry_selector(cell)(lattice))
    _assert_pattern_only(lattice, cell)
    if cell == "rand":
        assert len(_negative_edges(lattice)) == lattice.Ne_flips


@pytest.mark.parametrize("cell", ["rand", "randXERR"])
def test_erdos_renyi_nx_flips_once(cell, tmp_path):
    from kernels.ER import prepare_erdos_renyi

    args = SimpleNamespace(
        N=ER_N,
        p=ER_P,
        pflip=PFLIP,
        workdir=None,
        cell_type=cell,
        graph_engine="nx",
    )
    graph = prepare_erdos_renyi(args, path_data=tmp_path, path_plot=tmp_path)
    _assert_pattern_only(graph, cell)


def test_l3d_transcluster_lattice_flips_once():
    from kernels.L3D import make_transcluster_lattice

    args = SimpleNamespace(
        L=SIDE_3D, geometry="simple_cubic", p=PFLIP, workdir=None
    )
    lattice = make_transcluster_lattice(args, init_nw_dict=True)
    lattice.flip_sel_edges(lattice.nwDict["rand"]["G"])
    _assert_pattern_only(lattice, "rand")


@pytest.mark.parametrize("cell", ["rand", "randXERR"])
def test_load_or_compute_lattice2d_flips_once(cell, tmp_path):
    from lrgsglib.graphs.nx.Lattice2DNX.eigenspace import (
        load_or_compute_Lattice2DNX,
    )

    lattice = load_or_compute_Lattice2DNX(
        SIDE,
        "squared",
        cell_type=cell,
        save=False,
        compute="spectrum_numpy",
        pflip=PFLIP,
        init_nw_dict=True,
        path_data=tmp_path,
        path_plot=tmp_path,
    )
    _assert_pattern_only(lattice, cell)


# ---------------------------------------------------------------------------
# B2 - GT kernels draw an independent realization per call
# ---------------------------------------------------------------------------
def test_gt_lattice2d_realizations_are_independent():
    pytest.importorskip("graph_tool")
    from kernels.L2D import prepare_lattice

    args = SimpleNamespace(
        L=SIDE,
        geometry="squared",
        p=PFLIP,
        workdir=None,
        cell_type="rand",
        graph_engine="gt",
    )
    first = _edge_signature(prepare_lattice(args))
    second = _edge_signature(prepare_lattice(args))
    assert (first[2] < 0).any()
    assert not np.array_equal(first, second)


def test_gt_erdos_renyi_realizations_are_independent():
    pytest.importorskip("graph_tool")
    from kernels.ER import prepare_erdos_renyi

    args = SimpleNamespace(
        N=ER_N,
        p=ER_P,
        pflip=PFLIP,
        workdir=None,
        cell_type="rand",
        graph_engine="gt",
    )
    first = _edge_signature(prepare_erdos_renyi(args))
    second = _edge_signature(prepare_erdos_renyi(args))
    same = first.shape == second.shape and np.array_equal(first, second)
    assert not same


# ---------------------------------------------------------------------------
# B5 - concurrent sidecar writers
# ---------------------------------------------------------------------------
def test_concurrent_cfg_sidecar_writers_do_not_collide(tmp_path):
    from lrgsglib.statsys._naming import write_cfg_sidecar

    barrier = threading.Barrier(N_WRITERS)

    def write(i: int):
        barrier.wait()
        return write_cfg_sidecar(tmp_path, {"i": i}, {"model": 1}, {"run": 2})

    with ThreadPoolExecutor(N_WRITERS) as pool:
        futures = [pool.submit(write, i) for i in range(N_WRITERS)]
        errors = [f.exception() for f in futures]

    assert errors == [None] * N_WRITERS
    payload = json.loads((tmp_path / "cfg.json").read_text())
    assert set(payload) == {"graph", "model", "run"}
    assert [p.name for p in tmp_path.iterdir()] == ["cfg.json"]


# ---------------------------------------------------------------------------
# B7 - float formatting keeps grid points distinct
# ---------------------------------------------------------------------------
def test_run_dirname_floats_keep_grid_points_distinct():
    from lrgsglib.statsys._naming import format_value

    grid = np.linspace(0.1, 0.12, 41)  # %.3g collapsed it to 21 names
    assert len({format_value(float(p)) for p in grid}) == grid.size
    assert format_value(0.1 + 0.2) == "0.3"


def test_sbatch_array_script_passes_exact_p():
    from kernels.Serializer import render_sbatch_array_script

    script = render_sbatch_array_script(
        jobname="job",
        array_range=(0, 9),
        concurrent=2,
        mem_mb=1024,
        partition="",
        time="",
        exec_path="src/L3D_IsingDynamics.py",
        L=SIDE,
        p=0.1025,
        passthrough=[],
    )
    assert f" {SIDE} -p 0.1025 " in script


# ---------------------------------------------------------------------------
# B9 - temp names are independent of the stdlib random stream
# ---------------------------------------------------------------------------
def test_tmp_sibling_ignores_reseeded_random(tmp_path):
    from lrgsglib.statsys._io import tmp_sibling

    target = tmp_path / "cfg.json"
    random.seed(0)
    first = tmp_sibling(target)
    random.seed(0)
    second = tmp_sibling(target)
    assert first != second
    assert first.parent == second.parent == tmp_path
