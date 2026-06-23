"""Tests for Lattice2D animation tools (engine-agnostic; shared NX + GT).

``Lattice2D{NX,GT}.animate.states`` / ``.animate.largest_cluster`` and the
builders in ``lrgsglib/graphs/_shared/animation/lattice2d.py`` render
spin-configuration and largest-cluster animations (inline JS-HTML or a raw
matplotlib animation, with an optional gif/mp4 save). They were lifted from
``ipy/nbs/VM/VoterModel_signed_intro.ipynb`` and now live in the engine-neutral
``lrgsglib.graphs._shared.animation`` package (structural, bound on both
engines), so the same implementation backs both lattice engines. Every test
below runs against both ``engine='nx'`` and ``engine='gt'``.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from lrgsglib.graphs import Lattice2D
from lrgsglib.graphs._shared.animation.lattice2d import (
    make_lattice2d_animation,
    make_lattice2d_cluster_animation,
)
from lrgsglib.utils.statsys import edge_sign_arrays

SIDE = 4
GREY = (0.82, 0.82, 0.82)
POS = (0.84, 0.19, 0.15)

# Kernel tests inspect the returned image without rendering the FuncAnimation,
# which matplotlib flags on GC; that is expected here (rendering is covered by
# the high-level method tests via ``to_jshtml`` / ``save``).
pytestmark = pytest.mark.filterwarnings(
    "ignore:Animation was deleted without rendering anything"
)


@pytest.fixture(params=["nx", "gt"])
def lat(request):
    engine = request.param
    if engine == "gt":
        pytest.importorskip("graph_tool")
    kw = dict(side1=SIDE, geo="sqr", pbc=True, engine=engine, pflip=0.0, seed=1)
    if engine == "nx":
        kw["init_nw_dict"] = False  # NX-only construction flag
    return Lattice2D(**kw)


@pytest.fixture
def states(lat):
    rng = np.random.default_rng(0)
    return [rng.choice([-1, 1], size=lat.N).astype(np.int8) for _ in range(5)]


def _ragged(lat):
    idx, signs = [], []
    for nd in range(lat.N):
        js, sg = [], []
        for j, w in lat.get_neighbors_with_weights(nd):
            js.append(j)
            sg.append(-1 if w < 0 else 1)
        idx.append(np.asarray(js, np.int64))
        signs.append(np.asarray(sg, np.int8))
    return idx, signs


def _has_grey(rgb):
    return bool(np.any(np.all(np.isclose(rgb.reshape(-1, 3), GREY), axis=1)))


# --------------------------------------------------------------------------- #
# kernels (no IPython dependency)                                              #
# --------------------------------------------------------------------------- #
def test_states_kernel_returns_animation(lat, states):
    fig, ax = plt.subplots()
    res = make_lattice2d_animation(lat, fig, ax, states, vmin=-1, vmax=1)
    assert hasattr(res.animation, "to_jshtml")
    assert tuple(np.asarray(res.image.get_array()).shape) == tuple(lat.syshape)
    plt.close(fig)


def test_cluster_kernel_aligned_has_no_grey(lat):
    idx, signs = _ragged(lat)
    b = edge_sign_arrays(signs, "rawspin")
    fig, ax = plt.subplots()
    res = make_lattice2d_cluster_animation(
        lat, fig, ax, [np.ones(lat.N, np.int8)], idx, b, pos_color=POS)
    rgb = np.asarray(res.image.get_array())
    assert np.allclose(rgb.reshape(-1, 3)[0], POS)   # +1 colour everywhere
    assert not _has_grey(rgb)                          # nothing greyed out
    plt.close(fig)


def test_cluster_kernel_split_has_grey(lat):
    idx, signs = _ragged(lat)
    b = edge_sign_arrays(signs, "rawspin")
    s = np.ones(lat.N, np.int8)
    s[: lat.N // 2] = -1                                # two equal half-lattice domains
    fig, ax = plt.subplots()
    res = make_lattice2d_cluster_animation(lat, fig, ax, [s], idx, b)
    assert _has_grey(np.asarray(res.image.get_array()))  # non-largest domain is grey
    plt.close(fig)


def test_cluster_kernel_rejects_empty_frames(lat):
    idx, signs = _ragged(lat)
    b = edge_sign_arrays(signs, "rawspin")
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        make_lattice2d_cluster_animation(lat, fig, ax, [], idx, b)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# high-level lattice methods                                                   #
# --------------------------------------------------------------------------- #
def test_methods_return_inline_html(lat, states):
    pytest.importorskip("IPython")
    h1 = lat.animate.states(states, n_frames=3, fps=5)
    h2 = lat.animate.largest_cluster(states, cluster_mode="rawspin", n_frames=3, fps=5)
    assert hasattr(h1, "data") and len(h1.data) > 0
    assert hasattr(h2, "data") and len(h2.data) > 0


def test_animate_states_saves_gif(lat, states, tmp_path):
    pytest.importorskip("PIL")            # pillow writer for .gif
    pytest.importorskip("IPython")
    out = tmp_path / "anim.gif"           # explicit path -> used as given
    lat.animate.states(states, n_frames=3, fps=5, save=out)
    assert out.exists() and out.stat().st_size > 0


def test_inline_false_returns_matplotlib_animation(lat, states):
    # inline=False yields the raw matplotlib animation (pure matplotlib), not an
    # IPython HTML object -- so no IPython dependency is required.
    anim = lat.animate.states(states, n_frames=3, fps=5, inline=False)
    assert hasattr(anim, "to_jshtml") and hasattr(anim, "save")
    assert not hasattr(anim, "data")     # not an IPython.display.HTML
    plt.close("all")


def test_bare_save_name_resolves_under_plot_dir(tmp_path):
    # A bare filename is written under the lattice's own PLOT tree
    # (path_plot/<structure>/...); an explicit/absolute path is respected as
    # given. Tested on the pure resolver so no real animation is rendered.
    from types import SimpleNamespace
    from lrgsglib.graphs._shared.animation._common import resolve_plot_path

    # No path_data/sgpathn -> resolves directly under path_plot.
    lat_like = SimpleNamespace(path_plot=tmp_path, path_data=None, path_sgdata=None)
    assert resolve_plot_path(lat_like, "run.gif") == tmp_path / "run.gif"
    explicit = tmp_path / "sub" / "x.gif"
    assert resolve_plot_path(lat_like, explicit) == explicit
    assert resolve_plot_path(lat_like, None) is None

    # With a structure under path_data, a bare name nests path_plot/<structure>/.
    struct = SimpleNamespace(
        path_data=tmp_path, path_plot=tmp_path / "plot",
        path_sgdata=tmp_path / "l2d_squared",
    )
    assert resolve_plot_path(struct, "m.gif") == tmp_path / "plot" / "l2d_squared" / "m.gif"
