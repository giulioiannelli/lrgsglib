import numpy as np
import os
import sys


do_print = sys.argv[2]
slanzarv_OPT = sys.argv[3]
onlygraph_OPT = sys.argv[4]
L = int(sys.argv[1])
p_l = np.concatenate([np.linspace(0.075, 0.095, 3), 
                     np.linspace(0.095, 0.15, 7),
                     np.linspace(0.15, 0.45, 5)])
T_l = np.linspace(2/(L**2)**(1/2), 0.75, 40)

for p in p_l:
    for T in T_l:
        slanzarv_STR = f"slanzarv -m 2400 --nomail --jobname \"LRGSG_{L}_{T:.4g}_{p:.3g}\"" if slanzarv_OPT == "--slanzarv" else ""
        the_string = slanzarv_STR + " python src/lattc2dsq_IsingDynamics_mClusters.py "\
            + f"{L} {T:.4g} {p:.3g} -nA 500 {onlygraph_OPT}"
        if do_print == "--print":
            print(the_string)
        else:
            os.system(the_string)
        if onlygraph_OPT:
            break