"""Interactive-environment surface: ``from lrgsglib.notebooks import *``.

This package is the single import that labs (``ipy/labs/``) and report
notebooks (``ipy/nbs/``) use as their header. It is split by what a
notebook author reaches for, and every submodule is re-exported here in
full (no ``__all__`` anywhere, so star-imports chain through):

- :mod:`.session`        - root move, mpl style, paths, rng, verbosity shims,
                           filename helpers, library constants.
- :mod:`.common`         - kitchen-sink library chains (``shared``, ``core``,
                           ``plotlib``, ``utils``, ``utils.ipy``) plus raw
                           matplotlib artists, IPython display, ``json``,
                           ``floor``, ``py3Dmol``, ``plotly``.
- :mod:`.signed_graphs`  - engine-agnostic graph constructors and the cached
                           lattice loaders.
- :mod:`.spectral`       - LRG / spectral surface, reconstruction kernels,
                           protein-TMD reconstruction.
- :mod:`.dynamics`       - ``IsingDynamics``, SignedRW tables, CEM / SA
                           defaults and helpers, spin-configuration viewer.

The import order below is load-bearing: star-imports resolve name
collisions by "last one wins", and it mirrors the historical single-file
order. Keep it.

Canonical header (see ``.agents/rules/notebook-hygiene.md``)::

    from lrgsglib.notebooks import *
    use_lab_style()
    paths = make_lab_paths("<CODE>")
    seed  = None
    rng   = resolved_rng(seed)
"""
from .session import *
from .common import *
from .signed_graphs import *
from .spectral import *
from .dynamics import *

# Side effects (must stay last: `use_lab_style` needs the cwd that
# `move_to_rootf` establishes).
move_to_rootf()
use_lab_style()
