"""Base classes for the ``sg.plot`` / ``sg.animate`` accessors.

These give a signed graph a small, discoverable namespace for its visualisation
methods (``sg.plot.draw(...)``, ``sg.animate.voxel(...)``) instead of dozens of
methods on the class body. Per-graph-type accessors (e.g. ``Lattice2DAnimate``,
``Lattice3DAnimate``) subclass these and are exposed as ``plot`` / ``animate``
properties on the concrete graph classes.

The accessor is the single home for the figure-output policy: every method routes
its ``save`` through :func:`resolve_plot_path`, so animations/figures land under
``path_plot/<structure>/<subfolder>`` consistently (the subfolder derived from the
dynamics object that produced the data).
"""
from __future__ import annotations


def frames_and_model(source):
    """Split a positional ``source`` into ``(frames, model)``.

    A dynamics object (anything exposing an ``s_t`` trajectory -- duck-typed so
    this module never imports ``statsys``) contributes its trajectory as the
    frames *and* itself as the path-deriving model. Any other source (a sequence
    of state vectors / arrays) is treated as explicit frames with no model.
    """
    if hasattr(source, "s_t"):
        return source.s_t, source
    return source, None


class _Accessor:
    """Base for ``sg.plot`` / ``sg.animate``: holds a back-reference to the graph."""

    __slots__ = ("_sg",)

    def __init__(self, sg):
        self._sg = sg

    def __repr__(self):
        return f"{type(self).__name__}({type(self._sg).__name__})"
