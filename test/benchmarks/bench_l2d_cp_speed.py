import argparse
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from lrgsglib import ContactProcessEI, Lattice2D


def _make_log_steps(steps: int, nsamplelog: int) -> np.ndarray:
    """Return sorted unique sample steps including t=0."""
    if nsamplelog <= 0:
        return np.array([0], dtype=int)
    samples = np.geomspace(1, steps, nsamplelog, dtype=int)
    return np.unique(np.concatenate(([0], samples)))


def _run_py_with_log(cp: ContactProcessEI, steps: int, sample_steps: np.ndarray):
    """Run Python backend and record density at specified steps."""
    densities = np.empty(len(sample_steps), dtype=float)
    densities[0] = cp.s.mean()
    next_idx = 1

    start = perf_counter()
    for step in range(1, steps + 1):
        cp._sweep_ei_lambda_cache(  # type: ignore[attr-defined]
            cp.s,
            cp._lambda_arr,
            cp._neigh_indices,
            cp._neigh_weights,
            cp._neigh_offsets,
            cp._reverse_weights,
            cp.gamma_eff,
            cp._activation_fn,
            cp.N,
        )

        if next_idx < len(sample_steps) and step == sample_steps[next_idx]:
            densities[next_idx] = cp.s.mean()
            next_idx += 1

    elapsed = perf_counter() - start

    # If simulation ended before filling all slots, pad with final density
    while next_idx < len(sample_steps):
        densities[next_idx] = cp.s.mean()
        next_idx += 1

    return elapsed, densities


def run_once(
    L: int,
    pflip: float,
    gamma: float,
    steps: int,
    runlang: str,
    nsamplelog: int,
    clean_export: bool,
):
    """Run a single simulation and return timing, state snapshots, densities, and file info."""

    lattice = Lattice2D(L)
    lattice.flip_random_fract_edges(pflip)

    # Use nsamplelog for both backends (C uses num_log_samples).
    num_log_samples = nsamplelog if nsamplelog > 0 else 1000
    cp = ContactProcessEI(
        lattice,
        gamma=gamma,
        activation='relu',
        runlang=runlang,
        simref=1,
        ic='rand',
        num_log_samples=num_log_samples,
    )

    cp.init_contact_dynamics()
    s0 = cp.s_t[0].copy() if hasattr(cp, 's_t') and len(cp.s_t) > 0 else cp.s.copy()

    density_steps = None
    density_values = None
    dens_file = None

    if runlang.lower() == 'py' and nsamplelog > 0:
        sample_steps = _make_log_steps(steps, nsamplelog)
        elapsed, density_values = _run_py_with_log(cp, steps, sample_steps)
        density_steps = sample_steps.astype(float)
    else:
        start = perf_counter()
        cp.run(tqdm_on=False, steps=steps, verbose=False, clean_export=clean_export)
        elapsed = perf_counter() - start

        # For C backends, try to locate the exported density file
        if runlang.upper().startswith('C') and not clean_export:
            syshape = getattr(cp.sg, 'syshapePth', f"N={cp.sg.N}")
            dens_dir = Path(cp.sg.path_sgdata) / 'cntct' / syshape
            pattern = f"dens_*_{getattr(cp, 'out_id', '')}.bin"
            matches = sorted(dens_dir.glob(pattern))
            if matches:
                dens_file = matches[-1]
                density_values = np.fromfile(dens_file, dtype=np.float64)
                # Approximate the sampling steps for plotting
                density_steps = np.linspace(0, steps, num=len(density_values))

    s_final = cp.s.copy()

    return elapsed, s0, s_final, density_steps, density_values, dens_file


def main():
    parser = argparse.ArgumentParser(description="Benchmark ContactProcessEI Python vs C1c backends.")
    parser.add_argument('--L', type=int, default=64, help='Lattice size (default: 64)')
    parser.add_argument('--pflip', type=float, default=0.15, help='Edge flip fraction (default: 0.15)')
    parser.add_argument('--gamma', type=float, default=0.585, help='Infection rate gamma (default: 0.585)')
    parser.add_argument('--steps', type=int, default=100000, help='Simulation steps (default: 100000)')
    parser.add_argument('--runs', type=int, default=10, help='Number of runs to average timings (default: 10)')
    parser.add_argument('--nsamplelog', type=int, default=0,
                        help='Log-spaced samples for Python density logging (0 disables; also passed to C num_log_samples)')
    parser.add_argument('--clean-export', action='store_true',
                        help='Clean C-exported files (dens, snapshots). Default keeps them for plotting.')
    parser.add_argument('--plot', action='store_true', help='Generate comparison plot from the last run')
    args = parser.parse_args()

    print(f"Benchmarking with L={args.L}, pflip={args.pflip}, gamma={args.gamma}, steps={args.steps}, runs={args.runs}\n")
    if args.nsamplelog > 0:
        print(f"Python backend will log density at {args.nsamplelog} log-spaced samples.")
    if not args.clean_export:
        print("C backend exports will be kept (dens files recorded).")

    times_py = []
    times_c1c = []
    last_snapshots = None
    last_dens_py = None
    last_dens_c = None
    dens_files_c: list[Path] = []

    for i in range(args.runs):
        t_py, s0_py, s_final_py, steps_py, dens_py, _ = run_once(
            args.L, args.pflip, args.gamma, args.steps, 'py', args.nsamplelog, args.clean_export)
        t_c1c, s0_c1c, s_final_c1c, steps_c, dens_c, dens_path = run_once(
            args.L, args.pflip, args.gamma, args.steps, 'C1c', args.nsamplelog, args.clean_export)
        times_py.append(t_py)
        times_c1c.append(t_c1c)
        last_snapshots = (s0_py, s_final_py, s0_c1c, s_final_c1c)
        last_dens_py = (steps_py, dens_py)
        last_dens_c = (steps_c, dens_c)
        if dens_path:
            dens_files_c.append(dens_path)
        print(f"Run {i+1}/{args.runs}: Python={t_py:.3f}s, C1c={t_c1c:.3f}s, speedup={t_py/t_c1c if t_c1c else float('nan'):.2f}x")

    mean_py, std_py = np.mean(times_py), np.std(times_py)
    mean_c1c, std_c1c = np.mean(times_c1c), np.std(times_c1c)
    speedup = mean_py / mean_c1c if mean_c1c else float('nan')

    print("\nTiming summary (mean ± std over runs):")
    print(f"Python: {mean_py:.3f} ± {std_py:.3f} s")
    print(f"C1c:    {mean_c1c:.3f} ± {std_c1c:.3f} s")
    print(f"Speedup: {speedup:.2f}x")

    if dens_files_c:
        print("\nC density files saved:")
        for pth in dens_files_c:
            print(f"  {pth}")

    if args.plot and last_snapshots is not None:
        s0_py, s_final_py, s0_c1c, s_final_c1c = last_snapshots
        L = args.L
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))

        axes[0, 0].imshow(s0_py.reshape(L, L), cmap='binary', interpolation='nearest', vmin=0, vmax=1)
        axes[0, 0].set_title(f'Python: Initial (ρ={s0_py.mean():.3f})', fontsize=12, fontweight='bold')
        axes[0, 0].axis('off')

        axes[0, 1].imshow(s_final_py.reshape(L, L), cmap='binary', interpolation='nearest', vmin=0, vmax=1)
        axes[0, 1].set_title(f'Python: Final (ρ={s_final_py.mean():.3f})', fontsize=12, fontweight='bold')
        axes[0, 1].axis('off')

        axes[1, 0].imshow(s0_c1c.reshape(L, L), cmap='binary', interpolation='nearest', vmin=0, vmax=1)
        axes[1, 0].set_title(f'C1c: Initial (ρ={s0_c1c.mean():.3f})', fontsize=12, fontweight='bold')
        axes[1, 0].axis('off')

        axes[1, 1].imshow(s_final_c1c.reshape(L, L), cmap='binary', interpolation='nearest', vmin=0, vmax=1)
        axes[1, 1].set_title(f'C1c: Final (ρ={s_final_c1c.mean():.3f})', fontsize=12, fontweight='bold')
        axes[1, 1].axis('off')

        plt.suptitle(f'Contact Process: Python vs C1c\n(L={args.L}, pflip={args.pflip}, γ={args.gamma}, {args.steps} steps, ic=rand)',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        _fig_dir = Path(__file__).parent / ".tmp"
        _fig_dir.mkdir(exist_ok=True)
        plt.savefig(_fig_dir / 'contact_process_comparison.png', dpi=150, bbox_inches='tight')
        plt.close()
        print('\nPlot saved to /tmp/contact_process_comparison.png')

        print('\nState change analysis:')
        print(f'Python: {(s_final_py != s0_py).sum()}/{len(s0_py)} nodes changed ({100*(s_final_py != s0_py).sum()/len(s0_py):.1f}%)')
        print(f'C1c:    {(s_final_c1c != s0_c1c).sum()}/{len(s0_c1c)} nodes changed ({100*(s_final_c1c != s0_c1c).sum()/len(s0_c1c):.1f}%)')

        # Plot density vs time for the last run (Python backend)
        steps_py, dens_py = last_dens_py if last_dens_py else (None, None)
        if dens_py is not None and steps_py is not None:
            fig_py, ax_py = plt.subplots(figsize=(7, 4))
            ax_py.plot(steps_py, dens_py, 'o-', label='Python density')
            ax_py.set_xscale('log')
            ax_py.set_xlabel('Steps')
            ax_py.set_ylabel('Density')
            ax_py.legend()
            ax_py.set_title('Density vs time (Python backend)')
            plt.tight_layout()
            plt.close()
        elif args.nsamplelog > 0:
            print("Python density logging requested but no data captured.")

        # Plot density vs time for the last run (C1c backend)
        steps_c, dens_c = last_dens_c if last_dens_c else (None, None)
        if dens_c is not None and steps_c is not None:
            fig_c, ax_c = plt.subplots(figsize=(7, 4))
            ax_c.plot(steps_c, dens_c, 'o-', color='tab:red', label='C1c density')
            ax_c.set_xscale('log')
            ax_c.set_xlabel('Steps')
            ax_c.set_ylabel('Density')
            ax_c.legend()
            ax_c.set_title('Density vs time (C1c backend)')
            plt.tight_layout()
            plt.close()


if __name__ == "__main__":
    main()
