from lrgsglib.core import *

def main():
    # Setup
    move_to_rootf()
    plt.style.use('ipynb/nb_plotsheet.mplstyle')

    path_plotchladni = PATHPLOT / 'chladni'

    number_of_averages = 500
    side = 48
    geo = 'tri'
    pflip = 1
    kwargs_TriL2D = dict(
        side1=side, geo=geo, pflip=pflip, with_positions=True,
        path_data=path_plotchladni
    )

    # Lattice and basis setup
    ltmp = Lattice2D(**kwargs_TriL2D)
    ltmp.flip_random_fract_edges()
    ltmp.compute_k_eigvV(with_routine='numpy')

    basis = np.array([
        ltmp.get_eigV_check(i, reshaped=False)
        for i in range(ltmp.N // 2)
    ])

    # Dynamics parameters
    T_list = np.linspace(0., 2.2, 10)
    how_many = 3
    icc = 'uniform'
    runlang = 'C3'
    thrmSTEP = 1
    eqSTEP = 1
    freq = 5
    remove_files = True

    # Initialize storage
    spinovp = {T: [] for T in T_list}

    # Main loop
    for T in T_list:
        for _ in range(number_of_averages):
            kwargs_ISDY = dict(
                T=T, ic=icc, runlang=runlang,
                thrmSTEP=thrmSTEP, eqSTEP=eqSTEP,
                freq=freq, out_suffix=icc
            )
            isdy = IsingDynamics(ltmp, **kwargs_ISDY)
            isdy.init_ising_dynamics()
            isdy.run()

            if remove_files:
                isdy.remove_run_c_files(remove_stderr=True)
                ltmp.remove_exported_files()

            spinovp_tmp = np.array([
                spin_overlap(compute_recon(isdy.s, basis[:i + 1]), isdy.s)
                for i in range(ltmp.N // 2)
            ])
            spinovp[T].append(spinovp_tmp)
        save_path = ltmp.path_lrgsg / f'spinovp_{icc}_avg_std_{side}_{geo}_T{T:.3g}.pkl'
        save_path.parent.mkdir(parents=True, exist_ok=True)
        data = [
            np.mean(np.stack(spinovp[T]), axis=0),
            np.std(np.stack(spinovp[T]), axis=0)
        ]
        with open(save_path, 'wb') as f:
            pk.dump(data, f)
        print(f"Saved result to: {save_path}")

if __name__ == "__main__":
    main()
