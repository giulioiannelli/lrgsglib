"""Ising / SignedRW notebook surface: entry points, defaults, CEM & SA helpers.

The CEM / SA helpers compose `IsingDynamics` to produce standardized result
dicts and a 3D-slice spin visualizer. They are the building blocks for any
CHL/CEM-vs-SA report. No new physics primitives are introduced.
"""

from pathlib import Path as _Path

import numpy as _np

# IsingDynamics CEM / SA / topological defaults — surfaced at the
# notebook layer so labs and reports can configure realizations without
# importing deeper into the config tree.
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_CEM_ELITE_FRAC,
    DEFAULT_CEM_GREEDY,
    DEFAULT_CEM_GREEDY_SWEEPS,
    DEFAULT_CEM_INIT_SIGMA,
    DEFAULT_CEM_ITER,
    DEFAULT_CEM_POP_SIZE,
    DEFAULT_CEM_RESTARTS,
    DEFAULT_CEM_SIGMA_CEILING,
    DEFAULT_CEM_SIGMA_FLOOR,
    DEFAULT_CEM_SMOOTHING,
    DEFAULT_SA_COOLING_RATE,
    DEFAULT_SA_COOLING_SCHEDULE,
    DEFAULT_SA_N_TEMPERATURES,
    DEFAULT_SA_STEPS_PER_T,
)
from ..config.progargs.defs.IsingDynamics import DEFAULT_SA_STRONG_GREEDY_SWEEPS
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_GREEDY_SWEEPS as _SA_STRONG_GREEDY_SWEEPS,
)
from ..config.progargs.defs.IsingDynamics import DEFAULT_SA_STRONG_N_RESTARTS
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_N_RESTARTS as _SA_STRONG_N_RESTARTS,
)
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_N_TEMPERATURES,
)
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_N_TEMPERATURES as _SA_STRONG_N_TEMPERATURES,
)
from ..config.progargs.defs.IsingDynamics import DEFAULT_SA_STRONG_STEPS_PER_T
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_STEPS_PER_T as _SA_STRONG_STEPS_PER_T,
)
from ..config.progargs.defs.IsingDynamics import DEFAULT_SA_STRONG_T_FINAL
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_T_FINAL as _SA_STRONG_T_FINAL,
)
from ..config.progargs.defs.IsingDynamics import DEFAULT_SA_STRONG_T_INIT
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_STRONG_T_INIT as _SA_STRONG_T_INIT,
)
from ..config.progargs.defs.IsingDynamics import (
    DEFAULT_SA_T_FINAL,
    DEFAULT_SA_T_INIT,
    DEFAULT_TOPO_N_MODES,
    DEFAULT_TOPO_POLISH,
    DEFAULT_TOPO_POLISH_SWEEPS,
)

# IsingDynamics — direct interactive entry point for SA / CEM realizations.
from ..statsys.IsingDynamics import IsingDynamics

# SignedRW walker-side helpers (public API used by notebook diagnostics).
from ..statsys.SignedRW._kernel import signed_lattice_tables

# ======================================================================
# CEM hyperparameter exploration helpers (lab/notebook reuse)
# ======================================================================


def _strong_sa_default_settings() -> dict:
    """Default settings dict for `run_strong_sa`. Single source of truth."""
    return dict(
        T_init=_SA_STRONG_T_INIT,
        T_final=_SA_STRONG_T_FINAL,
        n_temperatures=_SA_STRONG_N_TEMPERATURES,
        steps_per_T=_SA_STRONG_STEPS_PER_T,
        n_restarts=_SA_STRONG_N_RESTARTS,
        greedy_sweeps=_SA_STRONG_GREEDY_SWEEPS,
    )


def _strong_sa_cache_key(sg, seed, settings) -> str:
    """Stable-hash cache key from graph signature, seed, and SA settings.

    The fingerprint includes a hash of the actual set of negative edges
    so that two `SignedGraph` instances with the same `(N, n_edges)` but
    different sign distributions (e.g. pre- and post-`flip_random_fract_edges`)
    do NOT collide on the same cache file.
    """
    import hashlib as _hashlib

    G = sg.gr["G"] if hasattr(sg, "gr") else None
    n_nodes = int(getattr(sg, "N", 0))
    n_edges = int(G.number_of_edges()) if G is not None else 0
    neg_sig = ""
    if G is not None:
        neg = sorted(
            ((u, v) if u <= v else (v, u))
            for u, v, d in G.edges(data=True)
            if d.get("weight", d.get("sign", 1.0)) < 0
        )
        neg_sig = _hashlib.sha1(repr(neg).encode()).hexdigest()[:12]
    payload = repr(
        (
            n_nodes,
            n_edges,
            neg_sig,
            int(seed) if seed is not None else None,
            sorted(settings.items()),
        )
    )
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


def run_strong_sa(
    sg,
    *,
    settings: dict | None = None,
    seed: int | None = None,
    cache_path: "_Path | None" = None,
    use_pybind: bool = True,
) -> dict:
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
            return {k: d[k] if d[k].ndim else d[k].item() for k in d.files} | {
                "settings": cfg,
                "seed": seed,
                "_cache_hit": True,
            }

    n_T = int(cfg["n_temperatures"])
    cooling_rate = (cfg["T_final"] / cfg["T_init"]) ** (1.0 / max(n_T - 1, 1))
    runlang = "pb_sa" if use_pybind else "py_sa"

    rng_master = _np.random.default_rng(seed)
    restart_seeds = rng_master.integers(
        0, 2**31 - 1, size=int(cfg["n_restarts"])
    )

    restart_E = _np.empty(int(cfg["n_restarts"]), dtype=_np.float64)
    restart_spins = _np.empty(
        (int(cfg["n_restarts"]), int(sg.N)), dtype=_np.int8
    )
    best_idx = -1
    best_E = _np.inf
    best_traj = None
    best_temps = None

    for r, sd in enumerate(restart_seeds):
        _np.random.seed(int(sd))  # BinDynSys.init_s uses np.random globally
        isdy = IsingDynamics(
            sg,
            ic="rand",
            runlang=runlang,
            sa_enabled=True,
            T_init=float(cfg["T_init"]),
            T_final=float(cfg["T_final"]),
            cooling_schedule="exponential",
            cooling_rate=float(cooling_rate),
            steps_per_T=int(cfg["steps_per_T"]),
            n_temperatures=int(cfg["n_temperatures"]),
            seed=int(sd),
        )
        isdy.init_ising_dynamics()
        isdy.run(sa_mode=True)
        _greedy_descent_at_T0(isdy, max_sweeps=int(cfg["greedy_sweeps"]))
        # Edge-sum energy per site, matching the C++ CEM kernel's
        # `persite_edge_energy` convention so SA and CEM are comparable.
        E_per_N = float(isdy.compute_energy()) / float(isdy.N)
        restart_E[r] = E_per_N
        restart_spins[r] = isdy.s.astype(_np.int8, copy=True)
        if E_per_N < best_E:
            best_E = E_per_N
            best_idx = r
            best_traj = (_np.asarray(isdy.sa_energy) / float(isdy.N)).copy()
            best_temps = _np.asarray(isdy.sa_temps).copy()

    out = {
        "final_E_per_N": best_E,
        "final_spins": restart_spins[best_idx].copy(),
        "restart_E": restart_E,
        "restart_spins": restart_spins,
        "traj_best": best_traj,
        "temps": best_temps,
        "settings": cfg,
        "seed": seed,
        "_cache_hit": False,
    }
    if cache_path is not None:
        _Path(cache_path).mkdir(parents=True, exist_ok=True)
        _np.savez(
            cache_file,
            final_E_per_N=_np.float64(best_E),
            final_spins=restart_spins[best_idx],
            restart_E=restart_E,
            restart_spins=restart_spins,
            traj_best=best_traj,
            temps=best_temps,
        )
    return out


def run_cem_pybind(
    sg, *, topo_n_modes: int, cem_kwargs: dict, seed: int | None = None
) -> dict:
    """Single pybind11 CEM realization on `sg`. Returns a standardized dict.

    `cem_kwargs` is forwarded as-is to `IsingDynamics.__init__`. Energies
    are returned per-site (the C++ kernel already divides by N).
    """
    isdy = IsingDynamics(
        sg,
        ic="rand",
        runlang="pb_topo_cem",
        topo_n_modes=int(topo_n_modes),
        seed=int(seed) if seed is not None else 0,
        **cem_kwargs,
    )
    isdy.init_ising_dynamics()
    isdy.run()
    return {
        "best_E_per_N": float(isdy.topo_cem_best_energy),
        "best_spins": _np.asarray(isdy.topo_cem_best_spins, dtype=_np.int8),
        "best_coeffs": _np.asarray(
            isdy.topo_cem_best_coeffs, dtype=_np.float64
        ),
        "restart_E": _np.asarray(
            isdy.topo_cem_restart_energies, dtype=_np.float64
        ),
        "history": _np.asarray(isdy.ene, dtype=_np.float64),
        "topo_n_modes": int(topo_n_modes),
        "cem_kwargs": dict(cem_kwargs),
        "seed": seed,
    }


def plot_lattice_slices_3d(
    spins,
    lattice_dim: tuple[int, int, int],
    *,
    axes=None,
    title: str | None = None,
    cmap: str = "bwr",
):
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
    axes[0].imshow(cube[:, :, mz], cmap=cmap, vmin=-1, vmax=1, origin="lower")
    axes[0].set_title(f"xy @ z={mz}")
    axes[1].imshow(cube[:, my, :], cmap=cmap, vmin=-1, vmax=1, origin="lower")
    axes[1].set_title(f"xz @ y={my}")
    axes[2].imshow(cube[mx, :, :], cmap=cmap, vmin=-1, vmax=1, origin="lower")
    axes[2].set_title(f"yz @ x={mx}")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    if title is not None:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes
