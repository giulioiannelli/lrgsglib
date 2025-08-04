# RBIsingM Programs

Code for simulating the random-bond Ising model.
Run `make c-make` from the repository root to compile all components
into `src/lrgsglib/Ccore/bin`.

Subdirectories:

- `base/` – C++ library exposing the `IsingModel` through Boost.Python.
- `simulatorC/` – standalone C executables (`IsingSimulator*.c`).
- `storer/` – utilities for serialising Ising configurations.
