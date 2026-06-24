"""Matplotlib 3D voxel-rendering primitives (engine- and graph-agnostic).

Pure plotting helpers for drawing a 3D boolean/colour array as a cube of solid
little voxels. ``Axes3D.voxels`` culls faces shared by two filled voxels, so a
packed cube would render as a hollow shell; the matplotlib "explode" trick
(upscale x2 with thin gaps) keeps every interior voxel a separate solid cube.
These know nothing about graphs -- they operate on plain numpy arrays. Used by
the structural 3D-lattice animations in
:mod:`lrgsglib.graphs._shared.animation.lattice3d`.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "VOX_KW",
    "VOX_ALPHA",
    "VOX_GAP",
    "explode",
    "gap_coords",
    "draw_voxels",
]

# Voxel cosmetics, shared across static views and movies.
VOX_KW = dict(edgecolor="w", lw=0.25)   # per-voxel edge borders
VOX_ALPHA = 0.5                          # default translucency of one phase
VOX_GAP = 0.05                           # gap so each cell is a separate cube


def explode(a):
    """Upscale x2 so neighbouring voxels no longer share (and thus cull) their
    common face: the originals land on even indices, the interleaved odd slots
    stay empty. Works on the ``(nx, ny, nz)`` mask and the ``(nx, ny, nz, 4)``
    colour array."""
    a = np.asarray(a)
    shp = np.array(a.shape)
    shp[:3] = shp[:3] * 2 - 1
    out = np.zeros(shp, dtype=a.dtype)
    out[::2, ::2, ::2] = a
    return out


def gap_coords(shape, *, gap=VOX_GAP):
    """Corner grid for an exploded voxel array, with thin ``gap`` gaps so every
    cell is drawn as a separate little cube instead of merging into a shell."""
    x, y, z = np.indices(np.array(shape) + 1).astype(float) // 2
    x[0::2, :, :] += gap
    x[1::2, :, :] += 1 - gap
    y[:, 0::2, :] += gap
    y[:, 1::2, :] += 1 - gap
    z[:, :, 0::2] += gap
    z[:, :, 1::2] += 1 - gap
    return x, y, z


def draw_voxels(ax, filled, rgba, *, gap=VOX_GAP, **kw):
    """Draw ``filled`` voxels coloured by ``rgba`` on a 3D axis, every cell a
    solid little cube (interior faces kept, not just the outer shell).

    ``filled`` is an ``(nx, ny, nz)`` boolean mask and ``rgba`` an
    ``(nx, ny, nz, 4)`` colour array (both un-exploded; exploded internally).
    Extra keyword args override :data:`VOX_KW` (edge cosmetics). Returns the
    ``ax.voxels`` collection dict (keyed by exploded voxel index), so callers
    can recolour faces in place for an animation.
    """
    fe, ce = explode(filled), explode(rgba)
    x, y, z = gap_coords(fe.shape, gap=gap)
    return ax.voxels(x, y, z, fe, facecolors=ce, **{**VOX_KW, **kw})
