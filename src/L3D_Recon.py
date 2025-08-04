# -*- coding: utf-8 -*-

from parsers.L3D_Recon import *
from kernels.L3D_Recon import *

def main():
    args = parse_arguments(parser)
    args_dict = vars(args)
    #
    locallog = initialize_custom_logger(args_dict, L3D_Recon_args)
    #
    compute_recon_prog_lattice_incr(args, locallog, prepare_lattice)

if __name__ == "__main__":
    main()


# if args.print_chrono:
#     Chronometer.enable()
# sovp, savepath = compute_reconstructions(args, locallog)
#
# save_results(sovp, args, savepath, locallog, prepare_lattice(args).syshapePth)
#
# if args.print_chrono:
#     print_accumulated_timings()