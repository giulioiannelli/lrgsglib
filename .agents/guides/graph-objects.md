# Graph Objects Guide

Quick reference for the graph object hierarchy, naming conventions, and
the `syshapePth` system used across the codebase.

## `syshapePth` — Graph Identity for File Paths

Every `SignedGraph` instance has a `syshapePth` string attribute that
uniquely identifies the graph's **size/shape** for use in directory
structures.  It is set at construction time and used to organize
simulation results.

### Format by Graph Type

| Graph Type | Format (typical) | Example | Determining params |
|------------|------------------|---------|-------------------|
| **Lattice2DNX** | `N={side1*side2}` (square) | `N=256` | `side1`, `side2` |
| | `L1={max}_L2={min}` (rect) | `L1=32_L2=16` | |
| | `...._prew={p}` (rewired) | `N=256_prew=0.1` | `prew` |
| **Lattice3DNX/GT** | `N={prod(dim)*mult}` (cubic) | `N=27000` | `dim`, `geo` |
| | `L0=_L1=_L2=_N=` (non-cubic) | `L0=30_L1=20_L2=10_N=6000` | |
| **ErdosRenyiNX/GT** | `N={n}_p={p:.3g}` | `N=500_p=0.05` | `n`, `p` |
| **MCG NX/GT** | `i={iter}_f={frac}_P={p1}_{p2}_{p3}_{p4}_{mode}` | `i=6_f=0.4_P=...` | cascade params |
| **FullyConnectedNX/GT** | `N={N}` | `N=100` | `N` |
| **WattsStrogatzNX** | `N={n}_p={p:.3g}` | `N=200_p=0.3` | `n`, `p` |
| **SierpinskiNX** | `N={syshape}_n={n}` | `N=3282_n=7` | `n` |
| **SignedGraph (base)** | `N={len(G)}` (fallback) | `N=42` | — |

### Where `syshapePth` is used

```
<path_sgdata>/
├── graph/<syshapePth>/       # Graph topology files
├── ising/<syshapePth>/       # Ising raw C output
├── ising_results/<syshapePth>/  # Checkpoint NPZ files  <-- ising kernels
├── spect/<syshapePth>/       # Spectral data
├── lrgsg/<syshapePth>/       # LRG flow data
└── phtra/<syshapePth>/       # Phase transition data
```

### NX vs GT implementation

- **NX**: Direct attribute `self.syshapePth = "N=..."` set in `__init__`
- **GT**: Private `self._syshapePth` with `@property` accessor + fallback

---

## Ising Result Filenames

Checkpoint files use the format:

```
<path_sgdata>/ising_results/<syshapePth>/<stem>_nt=<n_thermal>.npz
```

**Stem format**: `{prefix}_{runlang}_p={pflip:.3g}[_m={modes}]_q={idx:03d}`

**Examples:**
```
ising_results/N=27000/ene_pb_sa_p=0.25_q=001_nt=10.npz
ising_results/N=27000/magn_C3B_p=0.3_q=005_nt=10.npz
ising_results/L0=30_L1=20_L2=10_N=6000/ene_C1_p=0.1_q=001_nt=5.npz
```

The geometry (e.g. `simple_cubic`, `sqr`) is NOT in the filename — it's
already encoded in the parent path via `path_sgdata` (e.g. `l3d_simple_cubic/`).

The size/shape is NOT in the filename — it's the `syshapePth` subdirectory.

---

## Graph Constructor Patterns

### Lattice2D
```python
from lrgsglib.graphs.nx import Lattice2D
lat = Lattice2D(side1=64, geo='sqr', pflip=0.3)
lat.flip_random_fract_edges()
# lat.syshapePth → "N=4096"
```

### Lattice3D
```python
from lrgsglib.graphs.nx import Lattice3DNX
lat = Lattice3DNX(dim=30, geo='sc', pflip=0.2)
lat.flip_random_fract_edges()
# lat.syshapePth → "N=27000"
# lat.geo → "simple_cubic"  (normalized from 'sc')
```

### ErdosRenyi
```python
from lrgsglib.graphs.nx import ErdosRenyi
er = ErdosRenyi(n=500, p=0.05, pflip=0.1)
# er.syshapePth → "N=500_p=0.05"
```

---

## Key Design Rules

1. **`syshapePth` identifies size, not structure** — geometry/topology is
   in the parent directory name (via `path_sgdata`)
2. **Consistent across NX/GT** — same format for the same graph params
3. **Used for directories, not filenames** — filenames only contain
   run-specific params (runlang, pflip, quench, thermal count)
4. **Set once at construction** — do not modify after init
5. **Fallback**: base class defaults to `"N={num_nodes}"`

---

## File: `.agents/guides/graph-objects.md`
Last updated: 2026-03-13
