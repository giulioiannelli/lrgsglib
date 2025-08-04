from parsers.L2D_SlaplSpect import *

class EigenDistributionArgs:
    def __init__(self, workDir, navg, period, bins_count, mode, howmany, 
                 side, p, geo, eigmode=None, cell=None, verbose=False):
        self.workDir = workDir
        self.navg = navg
        self.period = period
        self.bins_count = bins_count
        self.mode = mode
        self.howmany = howmany
        self.side = side
        self.p = p
        self.geo = geo
        self.eigmode = eigmode
        self.cell = cell
        self.verbose = verbose

def process_eigen_distribution(fname_base, initial_data_fn, update_data_fn, 
                               save_data_fn, args):
    """
    Manages the overall process of computing and saving the eigenvalue or 
    eigenvector distributions. It handles reading existing data, updating 
    with new data, and saving results.

    Args:
        fname_base (str): Base name for the output file.
        initial_data_fn (function): Function to compute the initial data.
        update_data_fn (function): Function to update the bin counters.
        save_data_fn (function): Function to save the computed data.
        args (EigenDistributionArgs): Object containing all necessary params.
    """
    if args.verbose:
        print(f"Starting process_eigen_distribution for {fname_base}")

    # Determine the number of averages already done if a file exists
    navg_done = 0
    existing_files = sorted(glob.glob(f"{args.workDir}{fname_base}_*.pkl"))
    if existing_files:
        # Assuming the filename is of the form `fname_base_<navg_done>.pkl`
        navg_done = max(int(os.path.splitext(f.split('_')[-1])[0]) for f in existing_files)

    # If a file already exists, load the data to resume from where it left off
    if navg_done > 0:
        fname = args.workDir + fname_base + f"_{navg_done}.pkl"
        with open(fname, "rb") as f:
            bin_counter = pk.load(f)
        # Infer bins and bin_centers from existing data
        if args.mode == "eigvec_dist":
            initial_data = [val for counter in bin_counter for val in counter.elements()]
        else:
            initial_data = list(bin_counter.elements())
        if args.verbose:
            print(f"Loaded existing data from {fname}")
    else:
        fname = args.workDir + fname_base + f"_{navg_done}.pkl"
        initial_data = initial_data_fn(args)
        if len(initial_data) == 0:
            raise ValueError("Initial data is empty. Please check the input parameters or data generation.")
        if args.mode == "eigvec_dist":
            bin_counter = [Counter() for _ in range(args.howmany)]
        else:
            bin_counter = Counter()
    bins, bin_centers = create_symmetric_log_bins(initial_data, 
                                                args.bins_count)
    print(bins, bin_centers)
    nAvgNeed = args.navg - navg_done

    # Loop through periods until reaching the target number of averages
    for current_period in range((nAvgNeed // args.period) + bool(nAvgNeed % args.period)):
        if args.verbose:
            print(f"Processing period {current_period + 1}/{(nAvgNeed // args.period) + bool(nAvgNeed % args.period)}")
        
        # Calculate the batch size for this iteration
        batch_size = min(nAvgNeed - current_period * args.period, args.period)
        
        # Update the data using the batch size
        bin_counter = update_data_fn(
            batch_size, bins, bin_centers, bin_counter, args)
        
        # Save data with appropriate filename
        save_data_fn(args, bin_counter, fname)
        # Rename the file to reflect the updated number of averages
        navg_done += batch_size
        new_fname = args.workDir + fname_base + f"_{navg_done}.pkl"
        os.rename(fname, new_fname)
        fname = new_fname

    # Rename to the final state
    fname_final = args.workDir + fname_base + f"_{navg_done}.pkl"
    os.rename(fname, fname_final)
    if args.verbose:
        print(f"Renamed final file to {fname_final}")

def eigvec_initial_data(args):
    """
    Computes the initial eigenvector data for the lattice configuration.

    Args:
        args (EigenDistributionArgs): Object containing all necessary params.

    Returns:
        numpy.ndarray: The initial eigenvector values.
    """
    if args.verbose:
        print("Computing initial eigenvector data...")
    result = np.abs(eigV_for_lattice2D_ptch(side=args.side, pflip=args.p, 
                                          geo=args.geo, mode=args.eigmode, 
                                          howmany=args.howmany))
    if args.verbose:
        print(f"Computed initial eigenvector data of size {result.shape}")
    return result

def eigval_initial_data(args):
    """
    Computes the initial eigenvalue data for the lattice configuration.

    Args:
        args (EigenDistributionArgs): Object containing all necessary params.

    Returns:
        numpy.ndarray: The initial eigenvalue values.
    """
    if args.verbose:
        print("Computing initial eigenvalue data...")
    result = eigv_for_lattice2D(side=args.side, pflip=args.p, geo=args.geo)
    if args.verbose:
        print(f"Computed initial eigenvalue data of size {result.shape}")
    return result

def eigvec_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """
    Updates the bin counters for eigenvector values by adding data from 
    multiple samples.

    Args:
        batch_size (int): Number of samples to generate in this batch.
        bins (list): The bin edges.
        bin_centers (list): The bin centers.
        bin_counter (list): The current bin counters.
        args (EigenDistributionArgs): Object containing all necessary params.

    Returns:
        list: Updated bin_counter.
    """
    if args.verbose:
        print(f"Updating eigenvector data for batch size {batch_size}...")
    eig_values = [[] for _ in range(args.howmany)]
    for _ in range(batch_size):
        eigV = eigV_for_lattice2D_ptch(side=args.side, pflip=args.p, 
                                       geo=args.geo, mode=args.eigmode, 
                                       howmany=args.howmany)
        for i in range(args.howmany):
            eig_values[i].append(eigV[i])
    eig_values = [np.concatenate(i) for i in eig_values]

    if any(len(eig) == 0 for eig in eig_values):
        raise ValueError("Eigenvector data is empty. Please check the input parameters or data generation.")

    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = create_symmetric_log_bins(eig_values, args.bins_count)

    for i in range(args.howmany):
        bin_counter[i].update(bin_eigenvalues(eig_values[i], bins, 
                                              bin_centers))
    if args.verbose:
        print("Updated eigenvector data.")

    return bin_counter

def eigval_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """
    Updates the bin counters for eigenvalue values by adding data from 
    multiple samples.

    Args:
        batch_size (int): Number of samples to generate in this batch.
        bins (list): The bin edges.
        bin_centers (list): The bin centers.
        bin_counter (Counter): The current bin counter.
        args (EigenDistributionArgs): Object containing all necessary params.

    Returns:
        Counter: Updated bin_counter.
    """
    if args.verbose:
        print(f"Updating eigenvalue data for batch size {batch_size}...")
    eig_values = []
    for _ in range(batch_size):
        eigv = eigv_for_lattice2D(side=args.side, pflip=args.p, geo=args.geo)
        eig_values.append(eigv)
    eig_values = np.concatenate(eig_values)

    if len(eig_values) == 0:
        raise ValueError("Eigenvalue data is empty. Please check the input parameters or data generation.")

    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = create_symmetric_log_bins(eig_values, args.bins_count)

    bin_counter.update(bin_eigenvalues(eig_values, bins, bin_centers))
    if args.verbose:
        print("Updated eigenvalue data.")

    return bin_counter

def save_data(args, data, fname):
    """
    Saves the computed data to a file.

    Args:
        args (EigenDistributionArgs): Object containing all necessary params.
        data (list or Counter): Data to save.
        fname (str): The filename to save data into.
    """
    if args.verbose:
        print(f"Saving data to {fname}...")
    with open(fname, "wb") as f:
        pk.dump(data, f)
    if args.verbose:
        print(f"Data saved to {fname}")

def main():
    """
    Main function to parse arguments, create the necessary configurations, 
    and run the process for computing eigenvalue or eigenvector 
    distributions.
    """
    args = parser.parse_args()
    side = args.L
    p = args.p
    geo = args.geometry
    cell = args.cell_type
    mode = args.mode
    eigmode = args.eigen_mode
    navg = args.number_of_averages
    period = args.period
    bins_count = args.bins_count
    howmany = args.howmany
    verbose = args.verbose

    l = Lattice2D(side1=side, pflip=p, geo=geo, sgpathn=args.workDir)
    workDir = l.path_spect.name

    dist_args = EigenDistributionArgs(workDir, navg, period, bins_count, mode, 
                                      howmany, side, p, geo, eigmode, cell, verbose)

    match mode:
        case "eigvec_dist":
            fname_base = '_'.join([f"dist{howmany}", f"{p:.3g}", eigmode])
            process_eigen_distribution(fname_base, eigvec_initial_data,
                                       eigvec_update_data,
                                       save_data, dist_args)
        case 'eigval_dist':
            fname_base = '_'.join([f"dist_eigval", f"{p:.3g}", cell])
            process_eigen_distribution(fname_base, eigval_initial_data,
                                       eigval_update_data,
                                       save_data, dist_args)

if __name__ == "__main__":
    main()
