from parsers.ER_IsingDynamics import parser, parse_arguments
from kernels.ER_IsingDynamics import run_simulation
from kernels.ER import initialize_er_dict_args
from kernels.IsingDynamics import initialize_ising_dict_args, get_out_suffix
from lrgsglib import Chronometer


def main() -> None:
    """Main entry point for ER_IsingDynamics program."""
    args = parse_arguments(parser)

    # Determine if using ground_state initial condition
    ic_gs = args.init_cond.startswith('ground_state')
    number = int(args.init_cond.split('_')[-1]) if ic_gs else 0

    # Generate output suffix using merged logic from IsingDynamics
    out_suffix = get_out_suffix(args)

    # Initialize argument dictionaries
    erDictArgs = initialize_er_dict_args(args)
    isingDictArgs = initialize_ising_dict_args(args, out_suffix, args.NoClust)

    # Run simulation
    run_simulation(args, erDictArgs, isingDictArgs, number, args.remove_files)

    # Print timing information if requested
    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")


if __name__ == "__main__":
    main()
