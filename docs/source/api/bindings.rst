bindings - Language Bindings
============================

lrgsglib ships its performance-critical kernels as compiled C/C++ extensions
rather than as a single importable ``lrgsglib.bindings`` package. The native
sources are colocated with the statistical-physics models they accelerate, so
there is no standalone ``bindings`` Python module to autodocument here; the
code and its build machinery are described in the developer guide.

Where the native code lives
---------------------------

- **Shared C core** — reusable primitives (binary and continuous dynamical
  systems, clustering, continuous-time Markov chains, RNG helpers) live in
  ``lrgsglib/statsys/_ccore/`` (for example ``LRGSG_bindynsys.{c,h}``,
  ``LRGSG_clusters.{c,h}``, ``LRGSG_ctmc.{c,h}``).
- **Per-model kernels** — each model keeps its own compiled backend under
  ``lrgsglib/statsys/<Model>/ccore/``. These comprise plain C simulators (for
  example ``IsingMetropolis.c``, ``VoterSimulator.c``, ``ContactProcessEI.c``)
  and pybind11 modules that expose them to Python (for example
  ``ising_native.cpp``, ``voter_native.cpp``, and ``srw_native.cpp`` /
  ``random_walk.cpp`` under ``lrgsglib/statsys/SignedRW/ccore/``).

The compiled backends are selected at run time through the ``runlang=``
argument of each dynamics class (see :doc:`/api/statsys`), so user code never
imports the native modules directly.

See also
--------

- :doc:`/dev_guide/c_extensions` — writing and wrapping C/C++ kernels.
- :doc:`/dev_guide/build_system` — how the extensions are compiled.
- :doc:`/dev_guide/architecture` — where native code fits in the
  kernel → wrapper → API design.
