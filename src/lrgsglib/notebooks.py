from lrgsglib.config.funcs import move_to_rootf
from IPython.display import clear_output, display, HTML
import py3Dmol


from .shared import *
from .core import *
from .plotlib import *
from .utils.ipy import *

# Engine-agnostic graph factories (preferred interface in notebooks).
from .graphs import (
    SignedGraph,
    Lattice2D,
    Lattice3D,
    LatticeND,
    ErdosRenyi,
    StochasticBlockModel,
    HierarchicalModular,
    MultiplicativeCascade,
    VicsekGraph,
    DiracCombGraph,
    DiracBrushGraph,
    SierpinskiGraph,
    BarabasiAlbert,
    WattsStrogatz,
    FullyConnected,
    kRegularGraph,
    BipartiteGraph,
    RandomGeometric,
    ConfigurationModel,
    LFRBenchmark,
    HolmeKim,
    DualBarabasiAlbert,
    ExtendedBarabasiAlbert,
    DGMgraph,
)

# LRG / spectral utilities most commonly needed alongside SignedGraph methods.
from .utils.lrg import (
    get_graph_lspectrum,
    compute_entropy_observables_from_eigenvalues,
    specific_heat_tau_window,
    lapl_dists,
    extract_ultrametric_matrix,
    MakeLinkageMatrix,
    compute_normalized_linkage,
    compute_optimal_threshold,
    circular_layout_by_cluster,
    log_dendrogram,
    dendrogram_leaf_node_colors,
    compute_signed_diffusion_distance,
    compute_eigenmode_sign_distance,
    agmon_geodesic_distance,
)

# Spectral-reconstruction kernel (used by CHL Chladni-state notebooks).
from .utils.basic.linalg import compute_recon_ultra, compute_mse_from_recon

# Protein-TMD primitives (CHL-SF04 family — see utils/recon/protein/).
from .utils.recon.protein import (
    SubstrateSpec,
    build_substrate_basis,
    substrate_signature,
    spec_label,
    DEFAULT_SUBSTRATES,
    download_pdb,
    extract_ca_coordinates,
    extract_atoms_coordinates,
    coords_to_pdb_string_with_structure,
    assign_secondary_structure_from_coords,
    pad_protein_coordinates,
    safe_reconstruct_coordinates_from_features,
    protein_to_coordinate_feature_vector,
    kabsch_align,
    kabsch_rmsd,
    q3_score,
    # nb_helpers — used by the CHL-SF04 notebook
    load_protein_corpus,
    load_tmd_bundles,
    interp_mse_curve,
    mse_per_substrate,
    pflip_curves,
    auc_per_substrate,
    q3_curve_per_substrate,
    q3_emergence_per_substrate,
    ss_segments,
    layout_positions_2d,
    render_substrate_panel,
    reconstruct_protein,
)
import json   # surfaced to labs/notebooks via the canonical import-* pattern

# SignedRW walker-side helpers (public API used by notebook diagnostics).
from .statsys.SignedRW._kernel import signed_lattice_tables

# IsingDynamics — direct interactive entry point for SA / CEM realizations.
from .statsys.IsingDynamics import IsingDynamics

# Most-used library defaults (symbolic, reusable).
from .config.const import (
    PATHDATA,
    PATHPLOT,
    LRSG_ENTROPY_STEP,
    DEFAULT_ENTROPY_LEXPONENT,
    DEFAULT_ENTROPY_HEXPONENT,
    L2D_SIDE1,
    L2D_GEO_SQR,
    L2D_GEO_TRI,
    L2D_P_C_DICT,
    # CHL-SF04 protein-TMD experiment grid.
    PFLIP_SWEEP_GRID,
    K_FRAC_GALLERY,
    K_FRAC_Q3,
    N_SEEDS_PFLIP,
    CHL_SF04_MASTER_SEED,
    CHL_SF04_RES_MIN,
    CHL_SF04_RES_MAX,
    CHL_SF04_CORPUS_SIZE,
    CHL_SF04_N_FEATURES,
    CHL_SF04_LOG_MSE_FLOOR,
    Q3_EMERGENCE_THRESHOLD,
)

# IsingDynamics CEM / SA / topological defaults — surfaced at the
# notebook layer so labs and reports can configure realizations without
# importing deeper into the config tree.
from .config.progargs.defs.IsingDynamics import (
    DEFAULT_CEM_ITER,
    DEFAULT_CEM_POP_SIZE,
    DEFAULT_CEM_ELITE_FRAC,
    DEFAULT_CEM_INIT_SIGMA,
    DEFAULT_CEM_SMOOTHING,
    DEFAULT_CEM_SIGMA_FLOOR,
    DEFAULT_CEM_SIGMA_CEILING,
    DEFAULT_CEM_RESTARTS,
    DEFAULT_CEM_GREEDY,
    DEFAULT_CEM_GREEDY_SWEEPS,
    DEFAULT_TOPO_N_MODES,
    DEFAULT_TOPO_POLISH,
    DEFAULT_TOPO_POLISH_SWEEPS,
    DEFAULT_SA_T_INIT,
    DEFAULT_SA_T_FINAL,
    DEFAULT_SA_COOLING_SCHEDULE,
    DEFAULT_SA_COOLING_RATE,
    DEFAULT_SA_STEPS_PER_T,
    DEFAULT_SA_N_TEMPERATURES,
    DEFAULT_SA_STRONG_T_INIT,
    DEFAULT_SA_STRONG_T_FINAL,
    DEFAULT_SA_STRONG_N_TEMPERATURES,
    DEFAULT_SA_STRONG_STEPS_PER_T,
    DEFAULT_SA_STRONG_N_RESTARTS,
    DEFAULT_SA_STRONG_GREEDY_SWEEPS,
)

LINKAGE_METHOD = 'average'   # UPGMA — classical LRG choice
CMAP_CLUSTERS  = 'tab20'


# ======================================================================
# Phase 2 interactive-environment helpers (`from lrgsglib.notebooks import *`)
# ======================================================================
#
# These helpers back the canonical lab/report header described in
# `.agents/rules/lab-hygiene.md` and `.agents/rules/notebook-hygiene.md`:
#
#     from lrgsglib.notebooks import *
#     use_lab_style()
#     paths = make_lab_paths("<CODE>")
#     seed  = None
#     rng   = resolved_rng(seed)
#
# The `nprint` / `nlog` shims respect module-level flags so diagnostic
# output stays opt-in (default OFF).
#
# ----------------------------------------------------------------------

from pathlib import Path as _Path
import logging as _logging
import numpy as _np

verbose_print_nb: bool = False
verbose_log_nb: bool = False


def nprint(*args, **kwargs) -> None:
    """`print` shim: emits only when ``lrgsglib.notebooks.verbose_print_nb`` is True.

    Flip the flag in one cell to enable loud diagnostics within a lab
    or notebook session; leave it off everywhere else. Default OFF.
    """
    if verbose_print_nb:
        print(*args, **kwargs)


def nlog(msg, level: int = _logging.INFO, *args, **kwargs) -> None:
    """`logging` shim: emits only when ``lrgsglib.notebooks.verbose_log_nb`` is True.

    Logs under the ``lrgsg.nb`` logger. Default OFF.
    """
    if verbose_log_nb:
        _logging.getLogger("lrgsg.nb").log(level, msg, *args, **kwargs)


def resolved_rng(seed=None) -> "_np.random.Generator":
    """Return ``np.random.default_rng(seed)``. ``None`` → fresh entropy; int → reproducible.

    Use this exactly once per lab header so every stochastic call threads the
    same generator. Set ``seed`` to an int only when reproducing a specific
    result — hardcoding seeds in iterative exploration introduces systematic
    bias.
    """
    return _np.random.default_rng(seed)


def make_lab_paths(code_id: str) -> dict:
    """Canonical path dict for a project ``<code_id>``.

    Returns four sub-paths rooted at ``data/<code_id>/``:

        {"raw":     Path("data/<code_id>/raw"),
         "figures": Path("data/<code_id>/figures"),
         "cache":   Path("data/<code_id>/cache"),
         "cfg":     Path("data/<code_id>/cfg")}

    All parent directories are created on demand (``mkdir(parents=True,
    exist_ok=True)``). The returned dict is fresh on every call — callers
    may extend it with per-lab overrides.
    """
    root = _Path("data") / code_id
    out = {
        "raw":     root / "raw",
        "figures": root / "figures",
        "cache":   root / "cache",
        "cfg":     root / "cfg",
    }
    for p in out.values():
        p.mkdir(parents=True, exist_ok=True)
    return out


def use_lab_style() -> None:
    """Apply the interactive-environment mpl style sheet (``ipy/nb_plotsheet.mplstyle``).

    Safe to call multiple times; idempotent. Requires the current working
    directory to be the outer-repo root (``move_to_rootf`` handles this
    at module-import time).
    """
    import matplotlib.pyplot as _plt
    _plt.style.use("ipy/nb_plotsheet.mplstyle")


# ======================================================================
# CEM hyperparameter exploration helpers (lab/notebook reuse)
# ======================================================================
#
# These compose `IsingDynamics` to produce standardized result dicts and a
# 3D-slice spin visualizer. They are the building blocks for any future
# CHL/CEM-vs-SA report. No new physics primitives are introduced.

from .config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_T_INIT       as _SA_STRONG_T_INIT,
    DEFAULT_SA_STRONG_T_FINAL      as _SA_STRONG_T_FINAL,
    DEFAULT_SA_STRONG_N_TEMPERATURES as _SA_STRONG_N_TEMPERATURES,
    DEFAULT_SA_STRONG_STEPS_PER_T  as _SA_STRONG_STEPS_PER_T,
    DEFAULT_SA_STRONG_N_RESTARTS   as _SA_STRONG_N_RESTARTS,
    DEFAULT_SA_STRONG_GREEDY_SWEEPS as _SA_STRONG_GREEDY_SWEEPS,
)


def _strong_sa_default_settings() -> dict:
    """Default settings dict for `run_strong_sa`. Single source of truth."""
    return dict(
        T_init        = _SA_STRONG_T_INIT,
        T_final       = _SA_STRONG_T_FINAL,
        n_temperatures= _SA_STRONG_N_TEMPERATURES,
        steps_per_T   = _SA_STRONG_STEPS_PER_T,
        n_restarts    = _SA_STRONG_N_RESTARTS,
        greedy_sweeps = _SA_STRONG_GREEDY_SWEEPS,
    )


def _strong_sa_cache_key(sg, seed, settings) -> str:
    """Stable-hash cache key from graph signature, seed, and SA settings.

    The fingerprint includes a hash of the actual set of negative edges
    so that two `SignedGraph` instances with the same `(N, n_edges)` but
    different sign distributions (e.g. pre- and post-`flip_random_fract_edges`)
    do NOT collide on the same cache file.
    """
    import hashlib as _hashlib
    G = sg.gr['G'] if hasattr(sg, 'gr') else None
    n_nodes = int(getattr(sg, 'N', 0))
    n_edges = int(G.number_of_edges()) if G is not None else 0
    neg_sig = ''
    if G is not None:
        neg = sorted(((u, v) if u <= v else (v, u))
                     for u, v, d in G.edges(data=True)
                     if d.get('weight', d.get('sign', 1.0)) < 0)
        neg_sig = _hashlib.sha1(repr(neg).encode()).hexdigest()[:12]
    payload = repr((n_nodes, n_edges, neg_sig,
                    int(seed) if seed is not None else None,
                    sorted(settings.items())))
    return _hashlib.sha1(payload.encode()).hexdigest()[:16]


def _greedy_descent_at_T0(isdy, max_sweeps: int) -> int:
    """Deterministic single-flip descent: only flips spins with DeltaE < 0.

    Bypasses `IsingDynamics.metropolis` which calls `boltzmann_factor` and
    raises on T=0. Returns the number of sweeps actually performed (early
    exit at the first sweep with zero flips, i.e. local minimum).
    """
    sweeps_done = 0
    for _ in range(max_sweeps):
        flips = 0
        for node in range(isdy.N):
            neigh = isdy.neigh_wghtmagn(node)
            neigh_eng = isdy.neigh_ene(neigh)
            DeltaE = 2 * isdy.s[node] * neigh_eng
            if DeltaE < 0:
                isdy.flip_spin(node)
                flips += 1
        sweeps_done += 1
        if flips == 0:
            break
    return sweeps_done


def run_strong_sa(sg, *, settings: dict | None = None,
                  seed: int | None = None,
                  cache_path: "_Path | None" = None,
                  use_pybind: bool = True) -> dict:
    """Best-of-N strong simulated annealing on a SignedGraph.

    Geometric cooling from `T_init` to `T_final` over `n_temperatures`
    temperatures with `steps_per_T` Metropolis sweeps each, repeated for
    `n_restarts` independent runs from random init, each finished by a
    deterministic T=0 greedy descent of up to `greedy_sweeps` sweeps.
    Returns the best-of-N restart.

    Cached on disk under `cache_path / sa_strong_<key>.npz` when
    `cache_path` is given; key derives from graph signature + seed +
    settings hash.

    Returns
    -------
    dict
        ``{'final_E_per_N', 'final_spins', 'restart_E', 'restart_spins',
           'traj_best', 'temps', 'settings', 'seed'}``.
    """
    cfg = _strong_sa_default_settings()
    if settings:
        cfg.update(settings)

    if cache_path is not None:
        cache_key = _strong_sa_cache_key(sg, seed, cfg)
        cache_file = _Path(cache_path) / f"sa_strong_{cache_key}.npz"
        if cache_file.exists():
            d = _np.load(cache_file, allow_pickle=False)
            return {k: d[k] if d[k].ndim else d[k].item()
                    for k in d.files} | {'settings': cfg, 'seed': seed,
                                         '_cache_hit': True}

    n_T = int(cfg['n_temperatures'])
    cooling_rate = (cfg['T_final'] / cfg['T_init']) ** (1.0 / max(n_T - 1, 1))
    runlang = 'pb_sa' if use_pybind else 'py_sa'

    rng_master = _np.random.default_rng(seed)
    restart_seeds = rng_master.integers(0, 2**31 - 1,
                                        size=int(cfg['n_restarts']))

    restart_E      = _np.empty(int(cfg['n_restarts']), dtype=_np.float64)
    restart_spins  = _np.empty((int(cfg['n_restarts']), int(sg.N)),
                               dtype=_np.int8)
    best_idx       = -1
    best_E         = _np.inf
    best_traj      = None
    best_temps     = None

    for r, sd in enumerate(restart_seeds):
        _np.random.seed(int(sd))  # BinDynSys.init_s uses np.random globally
        isdy = IsingDynamics(
            sg,
            ic='rand',
            runlang=runlang,
            sa_enabled=True,
            T_init=float(cfg['T_init']),
            T_final=float(cfg['T_final']),
            cooling_schedule='exponential',
            cooling_rate=float(cooling_rate),
            steps_per_T=int(cfg['steps_per_T']),
            n_temperatures=int(cfg['n_temperatures']),
            seed=int(sd),
        )
        isdy.init_ising_dynamics()
        isdy.run(sa_mode=True)
        _greedy_descent_at_T0(isdy, max_sweeps=int(cfg['greedy_sweeps']))
        # Edge-sum energy per site, matching the C++ CEM kernel's
        # `persite_edge_energy` convention so SA and CEM are comparable.
        E_per_N = float(isdy.compute_energy()) / float(isdy.N)
        restart_E[r]     = E_per_N
        restart_spins[r] = isdy.s.astype(_np.int8, copy=True)
        if E_per_N < best_E:
            best_E = E_per_N
            best_idx = r
            best_traj = (_np.asarray(isdy.sa_energy)
                         / float(isdy.N)).copy()
            best_temps = _np.asarray(isdy.sa_temps).copy()

    out = {
        'final_E_per_N' : best_E,
        'final_spins'   : restart_spins[best_idx].copy(),
        'restart_E'     : restart_E,
        'restart_spins' : restart_spins,
        'traj_best'     : best_traj,
        'temps'         : best_temps,
        'settings'      : cfg,
        'seed'          : seed,
        '_cache_hit'    : False,
    }
    if cache_path is not None:
        _Path(cache_path).mkdir(parents=True, exist_ok=True)
        _np.savez(cache_file,
                  final_E_per_N=_np.float64(best_E),
                  final_spins  =restart_spins[best_idx],
                  restart_E    =restart_E,
                  restart_spins=restart_spins,
                  traj_best    =best_traj,
                  temps        =best_temps)
    return out


def run_cem_pybind(sg, *, topo_n_modes: int, cem_kwargs: dict,
                   seed: int | None = None) -> dict:
    """Single pybind11 CEM realization on `sg`. Returns a standardized dict.

    `cem_kwargs` is forwarded as-is to `IsingDynamics.__init__`. Energies
    are returned per-site (the C++ kernel already divides by N).
    """
    isdy = IsingDynamics(
        sg,
        ic='rand',
        runlang='pb_topo_cem',
        topo_n_modes=int(topo_n_modes),
        seed=int(seed) if seed is not None else 0,
        **cem_kwargs,
    )
    isdy.init_ising_dynamics()
    isdy.run()
    return {
        'best_E_per_N': float(isdy.topo_cem_best_energy),
        'best_spins'  : _np.asarray(isdy.topo_cem_best_spins, dtype=_np.int8),
        'best_coeffs' : _np.asarray(isdy.topo_cem_best_coeffs,
                                    dtype=_np.float64),
        'restart_E'   : _np.asarray(isdy.topo_cem_restart_energies,
                                    dtype=_np.float64),
        'history'     : _np.asarray(isdy.ene, dtype=_np.float64),
        'topo_n_modes': int(topo_n_modes),
        'cem_kwargs'  : dict(cem_kwargs),
        'seed'        : seed,
    }


def plot_lattice_slices_3d(spins, lattice_dim: tuple[int, int, int], *,
                           axes=None, title: str | None = None,
                           cmap: str = 'bwr'):
    """Render three orthogonal mid-slices of a 3D spin configuration.

    `spins` is a length-N array of {-1, +1} reshaped to `lattice_dim`
    (Lx, Ly, Lz). Three matplotlib `imshow` panels show the xy, xz, yz
    mid-planes. Returns `(fig, axes)`.
    """
    import matplotlib.pyplot as _plt
    Lx, Ly, Lz = lattice_dim
    cube = _np.asarray(spins).reshape(Lx, Ly, Lz)
    mx, my, mz = Lx // 2, Ly // 2, Lz // 2
    if axes is None:
        fig, axes = _plt.subplots(1, 3, figsize=(9, 3))
    else:
        fig = axes[0].figure
    axes[0].imshow(cube[:, :, mz], cmap=cmap, vmin=-1, vmax=1, origin='lower')
    axes[0].set_title(f'xy @ z={mz}')
    axes[1].imshow(cube[:, my, :], cmap=cmap, vmin=-1, vmax=1, origin='lower')
    axes[1].set_title(f'xz @ y={my}')
    axes[2].imshow(cube[mx, :, :], cmap=cmap, vmin=-1, vmax=1, origin='lower')
    axes[2].set_title(f'yz @ x={mx}')
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
    if title is not None:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes


move_to_rootf()
use_lab_style()
warnings.filterwarnings('ignore', category=scipy.sparse.SparseEfficiencyWarning)