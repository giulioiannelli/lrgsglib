import numpy as np
import os
import sys
#
noavg = 1000
List = [8, 16, 32, 64, 128, 256, 512]
plist = [0.01, 0.034, 0.08, 0.103, 0.12, 0.206, 0.412]
do_print = sys.argv[1]
slanzarv_OPT = sys.argv[2]
slanzarv_MEMMB = sys.argv[3]

for L in List:
    for p in plist:
        slanzarv_STR = (
            f'slanzarv -m {slanzarv_MEMMB} --nomail --jobname "lattc2dsq_perClusters_{L}_{p:.3g}"'
            if slanzarv_OPT == "--slanzarv"
            else ""
        )
        the_string = (
            slanzarv_STR
            + " python src/lattc2dsq_percClusters.py "
            + f"{L} {p:.3g} -nA {noavg} "
        )
        if do_print == "--print":
            print(the_string)
        else:
            print(the_string)
            os.system(the_string)
