"""Contact-process runner kernel: aggregation under the per-run layout.

kernels.ContactProcessDynamics._process_EI_C1c used to glob flat
``dens_*.bin`` files in the model directory — a pattern that matches
nothing under the Phase-C per-run-directory contract (the C binary now
writes ``rho_*`` inside the rundir and the model reads it back into the
density observable). These tests pin the fixed contract:

- the per-run series is collected off ``cp.density`` (no directory
  scan), so it works for every backend that populates the observable;
- the transient per-run directory is removed after collection;
- the cumulative ``dens_..._na=<i>.bin`` aggregate (kernel-owned naming,
  resume-compatible with existing data) accumulates across batches and
  keeps only the latest progress marker.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("networkx")

SRC = Path(__file__).resolve().parents[2] / "src"
if str(SRC) not in sys.path:  # kernels/ lives beside lrgsglib/ in src/
    sys.path.insert(0, str(SRC))

from kernels import ContactProcessDynamics as cpd


def _make_args(**over):
    base = dict(
        gamma=0.5,
        sp=None,
        init_cond="",
        cell_type="",
        out_suffix="",
        number_of_averages=2,
        save_frequency=1,
    )
    base.update(over)
    return SimpleNamespace(**base)


def _make_cp(pdir: Path, density, rundir_name: str):
    rundir = pdir / rundir_name
    rundir.mkdir(parents=True)
    (rundir / "cfg.json").write_text("{}")
    sg = SimpleNamespace(path_cntct=pdir, pflip=0.1)
    return SimpleNamespace(density=list(density), _c_rundir=rundir, sg=sg)


@pytest.fixture(autouse=True)
def _reset_kernel_state():
    cpd._batch_densities = []
    cpd._last_saved_index = 0
    yield
    cpd._batch_densities = []
    cpd._last_saved_index = 0


def test_process_collects_density_and_removes_rundir(tmp_path):
    args = _make_args()
    runs = [[0.5, 0.4, 0.3], [0.6, 0.5, 0.4]]
    for i, series in enumerate(runs, start=1):
        cp = _make_cp(tmp_path, series, f"p=0.1_gam=0.5_ns=3_lang=c_s={i}")
        args._current_average = i
        cpd._process_EI_C1c(cp, args)
        assert not cp._c_rundir.exists()  # transient per-run dir removed

    pdir, agg_prefix, sp_token = cpd._output_components(
        SimpleNamespace(path_cntct=tmp_path, pflip=0.1), args
    )
    # Only the latest progress marker survives.
    aggs = sorted(p.name for p in tmp_path.glob("dens_*_na=*.bin"))
    assert aggs == [f"dens_{agg_prefix}{sp_token}_na=2.bin"]
    data = np.fromfile(tmp_path / aggs[0], dtype=np.float64)
    assert np.allclose(data, np.concatenate([np.asarray(r) for r in runs]))
    # Resume scan finds the aggregate again.
    assert cpd._latest_saved_index(pdir, agg_prefix, sp_token) == 2


def test_process_skips_empty_density(tmp_path):
    """A run whose backend produced no density series contributes
    nothing (mirrors the old best-effort file hunt coming up empty)."""
    args = _make_args(number_of_averages=1)
    cp = _make_cp(tmp_path, [], "p=0.1_gam=0.5_ns=3_lang=c_s=1")
    args._current_average = 1
    cpd._process_EI_C1c(cp, args)
    assert not cp._c_rundir.exists()
    assert list(tmp_path.glob("dens_*.bin")) == []
