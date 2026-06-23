"""Shared helpers for the engine-agnostic graph plot/animation accessors."""
from __future__ import annotations

from pathlib import Path


def _relative_to(path, root):
    """``path`` relative to ``root``; ``None`` if either is missing or unrelated."""
    if path is None or root is None:
        return None
    try:
        return Path(path).relative_to(Path(root))
    except (ValueError, TypeError):
        return None


def resolve_plot_path(sg, save, *, model=None, subfolder=None):
    """Resolve a figure/animation ``save`` target under the graph's PLOT tree.

    Layout mirrors the data tree but rooted at ``path_plot``::

        path_plot / <structure> [/ <subfolder> [/ <size>]]

    where ``<structure>`` is the graph's own data subpath (``path_sgdata``
    relative to ``path_data``, e.g. ``l2d_squared``). The ``<subfolder>`` (and
    size, e.g. ``voter/N=...``) is *derived from a dynamics ``model``* -- its own
    ``dynpath`` re-rooted from the data tree onto the plot tree -- or given
    explicitly via ``subfolder``. With neither, figures land directly under
    ``path_plot / <structure>``.

    A bare filename lands under that directory; an absolute path, or one that
    already carries a directory component, is respected as-is. Returns ``None``
    when ``save`` is ``None``. Engine-robust: falls back to ``path_data/plot``
    when ``path_plot`` is unset (e.g. the GT engine).
    """
    if save is None:
        return None
    p = Path(save)
    if p.is_absolute() or p.parent != Path("."):
        return p  # explicit location -- respect it

    path_data = getattr(sg, "path_data", None)
    plot_root = getattr(sg, "path_plot", None)
    if plot_root is None and path_data is not None:
        plot_root = Path(path_data) / "plot"   # GT fallback (mirrors NX default)
    if plot_root is None:
        return p
    plot_root = Path(plot_root)

    rel = None
    if subfolder is not None:
        struct = _relative_to(getattr(sg, "path_sgdata", None), path_data)
        rel = (struct / str(subfolder)) if struct is not None else Path(str(subfolder))
    elif model is not None:
        # Mirror the dynamics' own data subtree (<structure>/<dyn>/<size>).
        rel = _relative_to(getattr(model, "dynpath", None), path_data)
    if rel is None:
        rel = _relative_to(getattr(sg, "path_sgdata", None), path_data)

    base = plot_root / rel if rel is not None else plot_root
    return base / p


# Back-compat: data-tree resolver (kept for any direct callers; the plot tree is
# the canonical home for figures/animations -- see ``resolve_plot_path``).
def resolve_save_path(lattice, save):
    """Resolve ``save`` against ``lattice.path_sgdata`` (the *data* tree)."""
    if save is None:
        return None
    p = Path(save)
    if p.is_absolute() or p.parent != Path("."):
        return p
    base = getattr(lattice, "path_sgdata", None)
    if base is None:
        return p
    return Path(base) / p
