# Lattice2D Spectral Analysis Test

## Overview

This test suite (`test_lattice2d_spectral.py`) demonstrates the complete functionality of the `Lattice2D` class with comprehensive spectral computations on signed graphs.

## What is Tested

### 1. **Basic Lattice Generation** ✓
- Creates lattices with different geometries:
  - **Squared** (coordination number z=4)
  - **Triangular** (coordination number z=6)
  - **Hexagonal** (coordination number z=3)
- Tests with periodic boundary conditions (PBC)
- Verifies signed edge generation with specified `pflip` (probability of negative edges)

### 2. **Spectral Computation** ✓
- Computes the full signed Laplacian spectrum
- Eigenvalue decomposition: `L_s = D_s - A`
  - `D_s`: Signed degree matrix (absolute degree values)
  - `A`: Signed adjacency matrix
- Verifies eigenvalue and eigenvector dimensions
- Analyzes ground state properties

### 3. **Eigenvector Analysis** ✓
- Retrieves individual eigenvectors
- Tests eigenvector binarization (converting continuous values to ±1)
- Analyzes node partitioning based on eigenvector signs
- Useful for identifying frustrated regions in spin-glass-like systems

### 4. **Clustering Analysis** ✓
- Performs connected component analysis from eigenvector partitioning
- Computes cluster size distributions
- Calculates the **infinite cluster probability** (P∞)
  - Percolation order parameter
  - Indicates phase transitions in frustrated systems

### 5. **Energy Computation** ✓
- **RBIM Energy** (Random Bond Ising Model)
  - Pairwise energy: `E = -∑_{<i,j>} w_{ij} s_i s_j`
  - Measures frustration in the system
- **Spherical SK Energy** (Sherrington-Kirkpatrick)
  - Mean-field spin-glass energy
  - Normalized by system size

### 6. **Edge Flipping** ✓
- Tests dynamic edge sign flipping
- Starts with unfrustrated system (pflip=0)
- Randomly flips edges to introduce frustration
- Verifies edge statistics after flipping

### 7. **Visualization** ✓
- Creates heatmap visualizations of eigenvectors
- Shows both continuous and binarized eigenvector values
- Saves output to `test/output/lattice2d_eigenvector_visualization.png`

## Running the Test

```bash
cd /path/to/lrgsglib
python test/test_lattice2d_spectral.py
```

## Requirements

The test requires the following packages (from `lrgsgenv.yml` or `pyproject.toml`):
- numpy
- scipy
- networkx
- matplotlib
- lmfit
- cupy (optional, for GPU acceleration)

## Output

The test produces:
- Console output with detailed test results
- Visualization PNG file in `test/output/`
- Summary of all tests passed/failed

## Physical Interpretation

The signed graph represents a frustrated system where:
- **Positive edges** (+1): Ferromagnetic interactions (spins prefer to align)
- **Negative edges** (-1): Antiferromagnetic interactions (spins prefer to anti-align)
- **pflip parameter**: Degree of frustration in the system

The **signed Laplacian eigenvectors** reveal:
- Ground state configuration minimizing frustration
- Cluster structure in the frustrated system
- Critical behavior near percolation transitions

## Example Results

```
Squared lattice (20×20, pflip=0.3):
  - Nodes: 400
  - Edges: 800
  - Negative edges: 240 (30%)
  - Ground state eigenvalue: ≈0
  - Full spectrum computed

Energy computations:
  - RBIM energy: -87.0
  - SK energy: -1.45
```

## Notes

- Some tests may show warnings when no clusters are found (uniform ground state)
- This is expected behavior when the ground state eigenvector has uniform sign
- Edge flipping may use a fallback method if the primary method encounters issues
- Visualization requires matplotlib backend support

## References

1. Kirkpatrick, S., & Sherrington, D. (1978). "Infinite-ranged models of spin-glasses." Physical Review B, 17(11), 4384.
2. Cucuringu, M. (2016). "Sync-Rank: Robust Ranking..." IEEE Trans. on Network Science and Engineering.
