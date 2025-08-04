# Ccore Programs

This directory contains the C and C++ sources for performance-critical parts of `lrgsglib`.
Most programs can be built by running `make c-make` from the repository root, which places
executables and Python extensions in `src/lrgsglib/Ccore/bin`.

- `LRGSG_utils.c` and related headers provide shared helper routines.
- `statsys/` contains statistical physics simulations such as the Ising and voter models.

Refer to the README in each subdirectory for details.
