import numpy as np
import os
import sys
#
L = int(sys.argv[1])
do_print = sys.argv[2]
slanzarv_OPT = sys.argv[3]
slanzarv_MEMMB = sys.argv[4]
p_l = np.linspace(0, 1, 50)
# p_l = np.concatenate(
#     [
#         np.linspace(0.075, 0.095, 3),
#         np.linspace(0.095, 0.15, 7),
#         np.linspace(0.15, 0.45, 5),
#     ]
# )
for p in p_l:
    slanzarv_STR = (
        f'slanzarv -m {slanzarv_MEMMB} --nomail --jobname "lattc2dsq_avgC_{L}_{p:.3g}"'
        if slanzarv_OPT == "--slanzarv"
        else ""
    )
    the_string = (
        slanzarv_STR
        + " python src/lattc2dsq_avgC.py "
        + f"{L} {p:.3g} -nA 256 "
    )
    if do_print == "--print":
        print(the_string)
    else:
        print(the_string)
        os.system(the_string)
