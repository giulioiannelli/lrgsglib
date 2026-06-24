"""Figure-saving helper (engine- and graph-agnostic).

One call for the recurring ``fig.savefig(..., bbox_inches='tight',
transparent=True)`` idiom, plus the dual-format (PDF + PNG) export and
parent-directory creation that labs/notebooks otherwise re-implement inline. All
cosmetic defaults are named constants in :mod:`lrgsglib.config.const`
(``SAVEFIG_*``), so they are tuned in one place rather than per call site.
"""
from __future__ import annotations

from pathlib import Path

from ..config.const import (
    SAVEFIG_BBOX,
    SAVEFIG_DPI,
    SAVEFIG_FORMATS,
    SAVEFIG_TRANSPARENT,
)

__all__ = ["save_fig"]


def save_fig(
    fig,
    path,
    *,
    formats=SAVEFIG_FORMATS,
    dpi=SAVEFIG_DPI,
    transparent=SAVEFIG_TRANSPARENT,
    bbox_inches=SAVEFIG_BBOX,
    mkdir=True,
    close=False,
):
    """Save ``fig`` to ``path``, defaulting to the lab's standard export style.

    ``path`` with an explicit suffix (``"fig.png"``) is written as-is, a single
    file. An *extensionless* ``path`` (``"fig"``) is written once per entry in
    ``formats`` -- by default ``fig.pdf`` and ``fig.png`` -- the common
    "vector for the paper, raster for the slide" pair. Parent directories are
    created when ``mkdir`` is True. Pass ``close=True`` to free the figure after
    writing. Returns the list of :class:`~pathlib.Path` objects written.

    All defaults (``dpi``, ``formats``, ``transparent``, ``bbox_inches``) come
    from the ``SAVEFIG_*`` constants in :mod:`lrgsglib.config.const`.
    """
    path = Path(path)
    if path.suffix:
        targets = [path]
    else:
        targets = [
            path.with_suffix(fmt if fmt.startswith(".") else f".{fmt}")
            for fmt in formats
        ]

    written = []
    for target in targets:
        if mkdir:
            target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            target, dpi=dpi, transparent=transparent, bbox_inches=bbox_inches
        )
        written.append(target)

    if close:
        import matplotlib.pyplot as plt

        plt.close(fig)

    return written
