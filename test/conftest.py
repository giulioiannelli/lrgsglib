"""
Root-level pytest configuration and fixtures for lrgsglib test suite.

Provides:
- ``--quick`` flag for fast smoke tests (small systems, minimal iterations)
- ``--show-plots`` flag for interactive plot display
- Parametrized fixtures for graph types, backends, and system sizes
- Helper for saving visual test outputs to test/.tmp/
"""

import shutil
from pathlib import Path

import matplotlib

import pytest
import numpy as np


# ---------------------------------------------------------------------------
# Collection exclusions
# ---------------------------------------------------------------------------

collect_ignore_glob = [
    "benchmarks/*",
    "scripts/*",
]


# ---------------------------------------------------------------------------
# CLI options
# ---------------------------------------------------------------------------

def pytest_addoption(parser):
    parser.addoption(
        "--quick",
        action="store_true",
        default=False,
        help="Use small systems for fast tests (<30s total)",
    )
    parser.addoption(
        "--show-plots",
        action="store_true",
        default=False,
        help="Display plots interactively (default: save only)",
    )
    parser.addoption(
        "--require-binaries",
        action="store_true",
        default=False,
        help=(
            "Fail (instead of skip) tests that need a compiled artifact "
            "(C binary, pybind/C++ extension). Use in built environments so "
            "a broken build cannot silently skip its test net."
        ),
    )


# ---------------------------------------------------------------------------
# Strict binary mode: convert missing-binary skips into failures
# ---------------------------------------------------------------------------

# Matches the skip reasons used by compiled-artifact guards across the suite
# ("<x> not built", "<x> binary not built", "native voter module unavailable",
# "C++ extension not available"). Optional *dependencies* (cupy, a GPU) keep
# skipping — only artifacts our own build produces are enforced.
import re as _re

_BINARY_SKIP_PATTERN = _re.compile(
    r"(?i)not built|binary|c\+\+ extension|native .*(unavailable|not available)"
)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if not (report.skipped and item.config.getoption("--require-binaries")):
        return
    if isinstance(report.longrepr, tuple):
        reason = str(report.longrepr[2])
    else:
        reason = str(report.longrepr)
    if _BINARY_SKIP_PATTERN.search(reason):
        report.outcome = "failed"
        report.longrepr = (
            f"--require-binaries: this test was about to SKIP because a "
            f"compiled artifact is missing ({reason}). Build it (make "
            f"cpp-make / the model's native target) or drop the flag."
        )


# ---------------------------------------------------------------------------
# Core mode fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def quick_mode(request):
    """Whether we're in quick mode."""
    return request.config.getoption("--quick")


@pytest.fixture
def show_plots(request):
    """Whether to display plots interactively."""
    return request.config.getoption("--show-plots")


@pytest.fixture
def system_size(quick_mode):
    """Lattice side length: 4 for quick, 16 for full."""
    return 4 if quick_mode else 16


@pytest.fixture
def n_steps(quick_mode):
    """Simulation steps: 50 for quick, 500 for full."""
    return 50 if quick_mode else 500


# ---------------------------------------------------------------------------
# Output directory for visual tests
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path(__file__).parent / ".tmp"


@pytest.fixture
def output_dir():
    """Path to test/.tmp/ for saving visual outputs."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    return OUTPUT_DIR


@pytest.fixture
def save_figure(show_plots):
    """Return a helper function that saves (and optionally shows) a figure."""

    def _save(fig, path, show=None):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=100, bbox_inches="tight")
        if show if show is not None else show_plots:
            import matplotlib.pyplot as plt
            plt.show()
        else:
            import matplotlib.pyplot as plt
            plt.close(fig)

    return _save


# ---------------------------------------------------------------------------
# Backend fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(params=["numpy", "scipy"])
def spectral_backend(request):
    """Parametrize over spectral backends."""
    return request.param


# ---------------------------------------------------------------------------
# Graph fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def small_lattice_sqr(tmp_path, system_size):
    """Square lattice, no frustration, PBC."""
    from lrgsglib.graphs.nx import Lattice2DNX
    return Lattice2DNX(
        side1=system_size, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )


@pytest.fixture
def frustrated_lattice(tmp_path, system_size):
    """Square lattice with pflip=0.3, edges flipped."""
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=system_size, geo="sqr", pbc=True, pflip=0.3,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    lat.flip_random_fract_edges()
    return lat


@pytest.fixture
def small_mcg(tmp_path):
    """Small MultiplicativeCascadeGraph for spectral tests."""
    from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX
    return MultiplicativeCascadeGraphNX(
        p1=0.8, p2=0.6, fraction=0.4, iterations=3,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )


@pytest.fixture
def small_er(tmp_path):
    """Small Erdos-Renyi graph."""
    from lrgsglib.graphs.nx import ErdosRenyiNX
    return ErdosRenyiNX(
        n=30, p=0.3, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
