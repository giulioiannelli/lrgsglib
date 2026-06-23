"""Engine-agnostic, graph-type-specific animation renderers.

These bind as methods on the concrete graph classes (mirroring
``graphs/_shared/_draw.py``): the *same* implementation drives both the
NetworkX and graph-tool engines, because it relies only on accessors present on
both (``syshape``, ``N``, ``get_neighbors_with_weights``). The animated quantity
is a generic state vector (e.g. the ``+/-1`` configuration of any ``BinDynSys``),
so the renderer depends on the graph *structure*, not on the dynamics that
produced the states. Output handling is delegated to the back-end-neutral
primitives in :mod:`lrgsglib.plotlib.animation`.

Graphs whose 2D rendering genuinely differs between engines bind their own
engine-specific renderers under ``graphs/<engine>/<Class>/animation.py`` instead.
"""
from .lattice2d import (
    animate_largest_cluster,
    animate_states,
    make_animation,
    make_lattice2d_animation,
    make_lattice2d_cluster_animation,
)
from .lattice3d import (
    animate_voxels,
    animate_voxels_cluster,
    show_voxels,
)

__all__ = [
    # 2D lattice (imshow movies)
    "make_animation",
    "animate_states",
    "animate_largest_cluster",
    "make_lattice2d_animation",
    "make_lattice2d_cluster_animation",
    # 3D lattice (voxel movies)
    "show_voxels",
    "animate_voxels",
    "animate_voxels_cluster",
]
