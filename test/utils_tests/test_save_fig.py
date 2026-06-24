"""Tests for the figure-saving helper :func:`lrgsglib.plotlib.save_fig`."""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from lrgsglib.config.const import PDF, PNG, SAVEFIG_FORMATS
from lrgsglib.plotlib import save_fig


@pytest.fixture
def fig():
    f, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    yield f
    plt.close(f)


def test_extensionless_writes_all_formats(fig, tmp_path):
    """A bare (extensionless) path writes one file per default format."""
    written = save_fig(fig, tmp_path / "myfig")
    assert {p.suffix for p in written} == set(SAVEFIG_FORMATS)
    for p in written:
        assert p.exists() and p.stat().st_size > 0
    assert (tmp_path / f"myfig{PDF}").exists()
    assert (tmp_path / f"myfig{PNG}").exists()


def test_explicit_extension_writes_single_file(fig, tmp_path):
    """An explicit suffix is honoured as-is: exactly one file."""
    written = save_fig(fig, tmp_path / "only.png")
    assert written == [tmp_path / "only.png"]
    assert written[0].exists()
    assert not (tmp_path / "only.pdf").exists()


def test_custom_formats(fig, tmp_path):
    written = save_fig(fig, tmp_path / "f", formats=(PNG,))
    assert written == [tmp_path / f"f{PNG}"]


def test_mkdir_creates_parents(fig, tmp_path):
    """Missing parent directories are created when mkdir=True."""
    target = tmp_path / "deep" / "nested" / "fig.png"
    written = save_fig(fig, target)
    assert written[0].exists()


def test_close_frees_figure(tmp_path):
    f, ax = plt.subplots()
    ax.plot([0, 1], [1, 0])
    num = f.number
    save_fig(f, tmp_path / "g.png", close=True)
    assert not plt.fignum_exists(num)
