from parsers.ER_IsingDynamics import parser,parse_arguments
from kernels.ER_IsingDynamics import *

def main():
    args = parse_arguments(parser)
    #
    ic_gs = args.init_cond.startswith('ground_state')
    number = int(args.init_cond.split('_')[-1]) if ic_gs else 0
    #
    out_suffix = get_out_suffix(args, ic_gs, number)
    erDictArgs = initialize_er_dict_args(args)
    isingDictArgs = initialize_ising_dict_args(args, out_suffix)
    #
    run_simulation(args, erDictArgs, isingDictArgs, number, args.remove_files)
    #
    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")

if __name__ == "__main__":
    main()
