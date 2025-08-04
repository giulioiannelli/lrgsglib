from lrgsglib import *
from parsers.L3D_Recon import L3D_Recon_progName, L3D_Recon_progNameShrt
from parsers.L3D_Recon_srun import parse_arguments, parser
#
def main():
    #
    args, unknown = parser.parse_known_args()
    progn = L3D_Recon_progName
    progn_shrt = L3D_Recon_progNameShrt
    #
    exec_bool, prnt_bool = args.exec, args.print
    side_list = args.dim_list
    pflp_list = args.pflip_linsp
    temp_list = args.Temp_linsp
    #
    # Memory function definition
    if args.slanzarv_minMB == args.slanzarv_maxMB:
        memoryfunc = lambda *_: args.slanzarv_minMB
    else:
        def memoryfunc(x):
            hl_side = [min(side_list), max(side_list)]
            hl_memy = [args.slanzarv_minMB, args.slanzarv_maxMB]
            return np.interp(x, hl_side, hl_memy).astype(int)
    #
    # Operate if either print or execute is requested
    if exec_bool or prnt_bool:
        count = count_exe = 0
        def operate(L, p, T):
            nonlocal count, count_exe
            progargs = [str(L), f"{p:.3g}", f"{T:.3g}"]
            opts = []
            opts.append("-m")
            opts.append(f"{memoryfunc(L)}")
            if args.nomail:
                opts.append("--nomail")
            if args.short:
                opts.append("--short")
            if args.moretime:
                opts.append(f"--time {args.moretime}")
            opts.append("--jobname")
            opts.append(f"{join_non_empty('_', progn_shrt, args.slanzarv_id, *progargs)}")
            execpath = LRGSG_SRC.relative_to(Path.cwd()) / f'{progn}.py'
            cmd = ['python', str(execpath)] + progargs + unknown
            slanz_cmd = ["slanzarv"] + opts + cmd
            #
            if prnt_bool:
                print(' '.join(slanz_cmd))
                count += 1
            if exec_bool:
                subprocess.run(slanz_cmd)
                count_exe += 1

        for L in side_list:
            for p in pflp_list:
                for T in temp_list:
                    operate(L, p, T)

        print(f"Total number of jobs executed: {count_exe}")
        print(f"Total number of jobs printed: {count}")

if __name__ == '__main__':
    main()