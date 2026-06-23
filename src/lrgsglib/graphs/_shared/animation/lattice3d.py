"""Engine-agnostic 3D-lattice voxel animations.

The 3D analogue of :mod:`lrgsglib.graphs._shared.animation.lattice2d`: where a 2D
configuration is shown with ``imshow(state.reshape(syshape))``, a 3D
configuration fills a cube of voxels coloured by the same length-``N`` state
vector reshaped to the 3-tuple ``lattice.syshape``. Like the 2D module it relies
only on accessors present on **both** engines (``syshape``, ``N``,
``get_neighbors_with_weights``), so ``Lattice3DNX`` and ``Lattice3DGT`` bind
these as methods (mirroring ``lattice2d.py``)::

    lat = Lattice3D(dim=12, geo='sc', engine='nx')   # or 'gt'
    lat.animate.voxel(states)                         # inline gif
    lat.animate.voxel(vm, save='run.gif')             # frames from vm.s_t; saved
                                                      # under <plot>/<structure>/voter/

The reusable voxel *geometry* primitives live in :mod:`lrgsglib.plotlib.voxels`
(plot-specific) and the rasterised-frames -> gif *codec* in
:mod:`lrgsglib.plotlib.animation.raster` (graph-agnostic); this module only owns
the lattice-aware orchestration. Only simple-cubic (``geo='sc'``) reshapes to a
clean cube.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ....plotlib.animation.raster import figure_to_image, frames_to_player_html
from ....plotlib.voxels import VOX_ALPHA, draw_voxels, explode
from ._accessor import _Accessor, frames_and_model
from ._common import resolve_plot_path

FrameLike = NDArray[np.integer] | NDArray[np.floating]

VOX_FRAMES = 60   # default number of frames kept for a voxel movie


def _state_cube(lattice, state) -> NDArray[np.floating]:
    """Reshape a length-``N`` state vector into the lattice's ``syshape`` cube."""
    syshape = getattr(lattice, "syshape", None)
    if syshape is None or len(tuple(syshape)) != 3:
        raise ValueError(
            "voxel animation needs a 3-tuple `syshape` on the lattice; got "
            f"{syshape!r}. Only simple-cubic Lattice3D exposes one (the GT "
            "engine derives it from `dim`)."
        )
    return np.asarray(state, dtype=float).reshape(*tuple(syshape))


def show_voxels(lattice, state, *, ax=None, cmap="hot", vmin=-1.0, vmax=1.0,
                alpha=VOX_ALPHA, alpha_sign=+1, title=None):
    """Static voxel view of one configuration (3D analogue of ``imshow``).

    The cube is drawn as solid little cubes (full volume, not just the shell);
    voxels whose value matches ``alpha_sign`` are translucent (``alpha``) so the
    interior stays visible. Pass ``alpha=1.0`` for a fully opaque cube.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    cube = _state_cube(lattice, state)
    cmo, norm = plt.get_cmap(cmap), Normalize(vmin, vmax)
    if ax is None:
        ax = plt.figure(figsize=(4, 4)).add_subplot(projection="3d")
    rgba = cmo(norm(cube))
    mid = 0.5 * (vmin + vmax)
    translucent = (cube >= mid) if alpha_sign > 0 else (cube < mid)
    rgba[..., 3] = np.where(translucent, alpha, 1.0)
    draw_voxels(ax, np.ones(cube.shape, bool), rgba)
    ax.set_axis_off()
    ax.set_aspect("equal")
    if title:
        ax.set_title(title)
    return ax


def _render_frames_local(rgba_cubes, syshape, recolor_edges, figsize, elev, azim):
    """Tessellate the full-cube geometry ONCE and recolour the per-voxel
    facecolours cube by cube (voxels can't blit), rasterising each to a PIL
    image. Shared verbatim by the serial path and by every parallel worker, so
    the frames are identical however they were produced. ``rgba_cubes`` is a
    sequence of ``(nx, ny, nz, 4)`` arrays; alpha-0 voxels are hidden (with
    ``recolor_edges`` their edges vanish too)."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection="3d")
    polys = draw_voxels(ax, np.ones(syshape, bool), rgba_cubes[0])
    ax.set_axis_off()
    ax.set_aspect("equal")
    if elev is not None or azim is not None:
        ax.view_init(elev=elev, azim=azim)
    for set_lim, n in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), syshape):
        set_lim(0, n)

    images = []
    for rgba in rgba_cubes:
        ce = explode(rgba)
        for key, coll in polys.items():
            fc = ce[key]
            coll.set_facecolor(fc)
            if recolor_edges:
                coll.set_edgecolor((1, 1, 1, 1.0 if fc[3] > 0 else 0.0))
        images.append(figure_to_image(fig))
    plt.close(fig)
    return images


def _render_frames_worker(rgba_cubes, syshape, recolor_edges, figsize, elev, azim):
    """Spawned-process entry point: pin a non-interactive backend (the worker is
    a fresh process) then render locally. Module-level so ``loky`` can pickle it
    by reference."""
    import matplotlib

    matplotlib.use("Agg", force=False)
    return _render_frames_local(rgba_cubes, syshape, recolor_edges, figsize, elev, azim)


def _resolve_n_jobs(n_jobs, n_frames):
    """Workers to use: ``None`` -> auto (<= 8 cores, only when there is enough
    work to amortise process startup), ``-1`` -> all cores, else the given count;
    falls back to serial (1) for short movies."""
    import os

    cpu = os.cpu_count() or 1
    workers = min(8, cpu) if n_jobs is None else (cpu if n_jobs < 0 else n_jobs)
    workers = max(1, min(workers, n_frames))
    return workers if (workers > 1 and n_frames >= 2 * workers) else 1


def _render_voxel_frames(rgba_cubes, syshape, *, recolor_edges, figsize, elev,
                         azim, n_jobs):
    """Render every coloured cube to an image, fanning the (independent,
    rasterisation-bound) frames across processes when it pays off and falling
    back to serial on any parallel failure. Round-robin chunks keep the load even
    and preserve frame order on reassembly."""
    workers = _resolve_n_jobs(n_jobs, len(rgba_cubes))
    if workers > 1:
        try:
            from joblib import Parallel, delayed

            chunks = [rgba_cubes[w::workers] for w in range(workers)]
            parts = Parallel(n_jobs=workers, backend="loky")(
                delayed(_render_frames_worker)(
                    c, syshape, recolor_edges, figsize, elev, azim
                )
                for c in chunks
            )
            images = [None] * len(rgba_cubes)
            for w, part in enumerate(parts):
                images[w::workers] = part
            return images
        except Exception:                       # spawn/joblib failure -> serial
            pass
    return _render_frames_local(rgba_cubes, syshape, recolor_edges, figsize, elev, azim)


def _voxel_movie(lattice, frames, color_of, *, fps, figsize, save,
                 recolor_edges=False, elev=None, azim=None, n_jobs=None):
    """Efficient voxel movie: the (cheap) per-frame colours are computed serially
    via ``color_of(s) -> (nx, ny, nz, 4)`` RGBA, then the (expensive) 3D frame
    rasterisation is fanned out across processes (``n_jobs``: ``None`` auto, ``1``
    serial, ``-1`` all cores) and the stack assembled into a looping gif (shown
    inline, optionally saved). Voxels with alpha 0 are hidden (set
    ``recolor_edges`` so their edges vanish too)."""
    syshape = tuple(getattr(lattice, "syshape"))
    rgba_cubes = [color_of(s) for s in frames]
    images = _render_voxel_frames(
        rgba_cubes, syshape, recolor_edges=recolor_edges, figsize=figsize,
        elev=elev, azim=azim, n_jobs=n_jobs,
    )
    return frames_to_player_html(images, fps=fps, save=save)


def animate_voxels(lattice, states, *, n_frames=VOX_FRAMES, fps=12, cmap="hot",
                   vmin=-1.0, vmax=1.0, figsize=(4, 4), save=None,
                   alpha=VOX_ALPHA, alpha_sign=+1, n_jobs=None, model=None,
                   subfolder=None):
    """3D voxel analogue of ``lat.animate_states`` (colour = spin).

    The cube is drawn as solid little cubes -- every interior voxel, not just
    the outer shell; voxels whose spin matches ``alpha_sign`` (default the +1
    phase) are rendered translucent with ``alpha`` while the opposite phase
    stays opaque and shows through. ``alpha`` is tunable; pass ``alpha=1.0`` for
    a fully opaque cube. ``n_frames`` evenly subsamples the trajectory; a bare
    ``save`` filename lands under ``lattice.path_sgdata``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    from ....utils.basic import subsample

    frames = subsample(list(states), n_frames) if n_frames else list(states)
    norm, cmo = Normalize(vmin, vmax), plt.get_cmap(cmap)
    mid = 0.5 * (vmin + vmax)

    def color_of(s):
        cube = _state_cube(lattice, s)
        rgba = cmo(norm(cube))
        translucent = (cube >= mid) if alpha_sign > 0 else (cube < mid)
        rgba[..., 3] = np.where(translucent, alpha, 1.0)
        return rgba

    return _voxel_movie(
        lattice, frames, color_of, fps=fps, figsize=figsize, n_jobs=n_jobs,
        save=resolve_plot_path(lattice, save, model=model, subfolder=subfolder),
    )


def animate_voxels_cluster(lattice, states, *, cluster_mode="rawspin",
                           n_frames=VOX_FRAMES, fps=12, figsize=(4, 4), save=None,
                           alpha=VOX_ALPHA, pos_color=(0.84, 0.19, 0.15),
                           neg_color=(0.13, 0.40, 0.74), n_jobs=None, model=None,
                           subfolder=None):
    """3D voxel analogue of ``lat.animate_largest_cluster``.

    Only the giant same-sign/satisfied domain is drawn, as solid little cubes
    coloured by its spin and rendered translucent (``alpha``, tunable; pass
    ``alpha=1.0`` => solid). The active-edge convention (``rawspin`` vs
    ``satisfied``) is read from this lattice's own signed adjacency.
    """
    from ....utils.basic import subsample
    from ....utils.statsys import cluster_components, signed_neighbor_arrays

    syshape = tuple(getattr(lattice, "syshape"))
    idx, b = signed_neighbor_arrays(lattice, cluster_mode)
    frames = subsample(list(states), n_frames) if n_frames else list(states)

    def color_of(s):
        s = np.asarray(s, np.int8).ravel()
        label = cluster_components(s, idx, b)
        largest = label == np.bincount(label).argmax()
        col = pos_color if (largest.any() and int(s[largest][0]) > 0) else neg_color
        rgba = np.zeros((*syshape, 4))                      # alpha 0 => hidden
        rgba[largest.reshape(*syshape)] = (*col, alpha)
        return rgba

    return _voxel_movie(
        lattice, frames, color_of, fps=fps, figsize=figsize, n_jobs=n_jobs,
        save=resolve_plot_path(lattice, save, model=model, subfolder=subfolder),
        recolor_edges=True,
    )


# === Accessors (bound as ``lat.plot`` / ``lat.animate`` on both engines) ===


class _Lattice3DAnimate(_Accessor):
    """``lat.animate`` for 3D lattices: voxel movies.

    The positional ``source`` may be a dynamics object (e.g. a ``VoterModel`` --
    its ``s_t`` trajectory is the frames and the gif is filed under that model's
    plot subtree) or an explicit sequence of state vectors.
    """

    def voxel(self, source, *, model=None, subfolder=None, **kw):
        frames, src_model = frames_and_model(source)
        return animate_voxels(
            self._sg, frames, model=model or src_model, subfolder=subfolder, **kw
        )

    def voxel_cluster(self, source, *, model=None, subfolder=None, **kw):
        frames, src_model = frames_and_model(source)
        return animate_voxels_cluster(
            self._sg, frames, model=model or src_model, subfolder=subfolder, **kw
        )


class _Lattice3DPlot(_Accessor):
    """``lat.plot`` for 3D lattices: static renderings."""

    def voxel(self, state, **kw):
        """Static voxel view of one configuration (3D analogue of ``imshow``)."""
        return show_voxels(self._sg, state, **kw)
