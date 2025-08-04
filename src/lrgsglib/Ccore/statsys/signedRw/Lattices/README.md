# Random Walker on Signed Lattices

`random_walker.cpp` implements a walker that tracks the sign accumulated
while moving on a sparse lattice. The module is built with pybind11.

Compile the extension using `make c-make`. The resulting shared object
is stored in `src/lrgsglib/Ccore/bin/`.
