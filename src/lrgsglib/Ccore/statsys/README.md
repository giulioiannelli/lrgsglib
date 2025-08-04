# Statistical Physics Programs

This directory hosts C and C++ implementations of various dynamical models.
The root `Makefile` provides a target `c-make` that builds all binaries and
extensions into `src/lrgsglib/Ccore/bin`.

Contents:

- `RBIsingM/` – programs for the random-bond Ising model.
- `signedRw/` – a random walker on signed lattices.
- `voterM/` – an implementation of the voter model.
- `contactP/` – a minimal contact-process simulator.
- `random_walk.cpp` – minimal pybind11 example of a random walk.

See each subdirectory for specific build notes.
