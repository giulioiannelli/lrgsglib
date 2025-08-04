from lrgsglib import *
from .L2D import *

def perform_spectral_calculations(args):
    """Main function to parse arguments and run the process."""
    # Determine the filename base based on the mode
    if args.mode.endswith("dist"):
        if args.mode == "eigvec_dist":
            fname_base = f"dist{args.howmany}_{args.p:.3g}_{args.eigen_mode}"
            initial_fn = eigvec_initial_data
            update_fn = eigvec_update_data
        elif args.mode == "eigval_dist":
            fname_base = f"dist_eigval_{args.p:.3g}_{args.cell_type}"
            initial_fn = eigval_initial_data
            update_fn = eigval_update_data
        # Process the eigen distribution
        process_eigen_distribution(fname_base, initial_fn, update_fn, save_data,
                                args)
    if args.mode == "eigvals":
        fname_base = f"eigvals_{args.p:.3g}_{args.cell_type}"
        #
        eigvlist = []
        for _ in range(args.number_of_averages):
            eigv = eigv_for_lattice2D(side=args.L, pflip=args.p, geo=args.geo, 
                           mode='_'.join(["some", str(args.howmany)]))
            eigvlist.append(eigv)
            if _ % args.period == 0:
                path_fname_base = Lattice2D(
                    side1=args.L, pflip=args.p, geo=args.geo, path_data=args.workDir
                ).path_spect
                path_fname = path_fname_base / Path(f"{fname_base}_{_}.pkl")
                path_fname_1 = path_fname_base / Path(f"{fname_base}_{_ - args.period}.pkl")
                if os.path.exists(path_fname_1):
                    os.remove(path_fname_1)
                with open(path_fname, "wb") as f:
                    pk.dump(eigvlist, f)

def process_eigen_distribution(
        fname_base, initial_data_fn, update_data_fn, save_data_fn, args):
    """
    Manages the overall process of computing and saving the eigenvalue or
    eigenvector distributions.
    """
    if args.verbose:
        print(f"Starting process_eigen_distribution for {fname_base}")

    # Check for existing files and determine the number of averages already done
    navg_done = 0
    working_path = Lattice2D(
        side1=args.L, pflip=args.p, geo=args.geo, path_data=args.workDir
    ).path_spect
    search_str = f"{fname_base}_*.pkl"
    existing_files = sorted(glob.glob(str(working_path / Path(search_str))))
    #
    if existing_files:
        navg_done = max(int(os.path.splitext(f.split('_')[-1])[0])
                        for f in existing_files)

    # Load existing data or compute initial data if no files exist
    if navg_done > 0:
        path_fname = working_path / Path(f"{fname_base}_{navg_done}.pkl")
        with open(path_fname, "rb") as f:
            bin_counter = pk.load(f)
        if args.mode == "eigvec_dist":
            initial_data = [val for counter in bin_counter for val in counter.elements()]
        else:
            initial_data = list(bin_counter.elements())
        if args.verbose:
            print(f"Loaded existing data from {path_fname}")
    else:
        path_fname = working_path / Path(f"{fname_base}_{navg_done}.pkl")
        initial_data = initial_data_fn(args)
        if initial_data is None:
            raise ValueError("Initial data is empty. Please check the input parameters or data generation.")
        if args.mode == "eigvec_dist":
            bin_counter = [Counter() for _ in range(args.howmany)]
        else:
            bin_counter = Counter()

    # Create bins based on mode
    if args.mode == "eigvec_dist":
        bins, bin_centers = create_symmetric_log_bins(initial_data, args.bins_count)
    else:
        bins, bin_centers = linear_binning_hist(initial_data, args.bins_count)
    nAvgNeed = args.number_of_averages - navg_done
    total_periods = (nAvgNeed // args.period) + bool(nAvgNeed % args.period)

    # Process each period until reaching the required number of averages
    for current_period in range(total_periods):
        if args.verbose:
            print(f"Processing period {current_period + 1}/{total_periods}")
        batch_size = min(nAvgNeed - current_period * args.period, args.period)
        bin_counter = update_data_fn(batch_size, bins, bin_centers,
                                     bin_counter, args)
        save_data_fn(args, bin_counter, path_fname)
        navg_done += batch_size
        new_fname = working_path /  Path(f"{fname_base}_{navg_done}.pkl")
        os.rename(path_fname, new_fname)
        path_fname = new_fname
    # Rename the final file
    fname_final = working_path /  Path(f"{fname_base}_{navg_done}.pkl")
    os.rename(path_fname, fname_final)
    if args.verbose:
        print(f"Renamed final file to {fname_final}")

def eigvec_initial_data(args):
    """Computes the initial eigenvector data."""
    if args.verbose:
        print("Computing initial eigenvector data...")
    result = np.abs(eigV_for_lattice2D_ptch(
        side=args.L, pflip=args.p, geo=args.geo, mode=args.eigen_mode,
        howmany=args.howmany))
    if args.verbose:
        print(f"Computed initial eigenvector data of size {result.shape}")
    return result

def eigval_initial_data(args):
    """Computes the initial eigenvalue data."""
    if args.verbose:
        print("Computing initial eigenvalue data...")
    result = eigv_for_lattice2D(
        side=args.L, pflip=args.p, geo=args.geo)
    if args.verbose:
        print(f"Computed initial eigenvalue data of size {result.shape}")
    return result

def eigvec_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """Updates the bin counters for eigenvector values."""
    if args.verbose:
        print(f"Updating eigenvector data for batch size {batch_size}...")
    eig_values = [[] for _ in range(args.howmany)]
    for _ in range(batch_size):
        eigV = eigV_for_lattice2D_ptch(
            side=args.L, pflip=args.p, geo=args.geo, mode=args.eigen_mode,
            howmany=args.howmany)
        for i in range(args.howmany):
            eig_values[i].append(eigV[i])
    eig_values = [np.concatenate(i) for i in eig_values]

    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = create_symmetric_log_bins(eig_values, args.bins_count)

    for i in range(args.howmany):
        bin_counter[i].update(bin_eigenvalues(eig_values[i], bins, bin_centers))
    if args.verbose:
        print("Updated eigenvector data.")
    return bin_counter

def eigval_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """Updates the bin counters for eigenvalue values."""
    if args.verbose:
        print(f"Updating eigenvalue data for batch size {batch_size}...")
    eig_values = [eigv_for_lattice2D(
        side=args.L, pflip=args.p, geo=args.geo) for _ in range(batch_size)]
    eig_values = np.concatenate(eig_values)

    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = linear_binning_hist(eig_values, args.bins_count)

    bin_counter.update(bin_eigenvalues(eig_values, bins, bin_centers))
    if args.verbose:
        print("Updated eigenvalue data.")
    return bin_counter

def compute_eigval(args):
    """Computes the eigenvalues for the lattice."""
    if args.verbose:
        print("Computing eigenvalues...")
    result = eigv_for_lattice2D(
        side=args.L, pflip=args.p, geo=args.geo)
    if args.verbose:
        print(f"Computed eigenvalues of size {result.shape}")
    return result

def save_data(args, data, fname):
    """Saves the computed data to a file."""
    if args.verbose:
        print(f"Saving data to {fname}...")
    with open(fname, "wb") as f:
        pk.dump(data, f)
    if args.verbose:
        print(f"Data saved to {fname}")