from parsers.L2D_SlaplSpect_Serialiser import *
#
args = parser.parse_args()
#
geo = args.geometry
cell = args.cell_type
mode = args.mode
eigmode = args.eigen_mode
navg = args.number_of_averages
T = args.period
bins_count = args.bins_count
howmany = args.howmany
#
progName = L2D_SlaplSpect_progName
progNameShrt = L2D_SlaplSpect_progNameShrt
progMode = '_'.join(args.mode.split('_')[1:])
execBool = args.exec
printBool = args.print
#
if mode.endswith("eigvec_dist"):
    List = 2 ** np.arange(4, 10)
    plist = np.linspace(0.06, 0.115, num=10)
if mode.endswith("eigval_dist"):
    List = [16, 32, 48, 64, 96, 128]
    plist = [0.01, 0.025, 0.05, 0.075, 0.08, 0.09, 0.1, 0.11, 0.125, 0.15, 0.25, 0.5, 0.75]
if mode.startswith("slanzarv"):
    if args.slanzarv_minMB == args.slanzarv_maxMB:

        def memoryfunc(*_):
            return args.slanzarv_minMB

    else:

        def memoryfunc(x):
            return int(
                np.interp(
                    x,
                    [min(List), max(List)],
                    [args.slanzarv_minMB, args.slanzarv_maxMB],
                )
            )

    def slanzarv_str(mode, L, p, geo, c):
        slanzarvopt = "--nomail --jobname "
        slanzarvstr = f"slanzarv -m {memoryfunc(L)} {slanzarvopt}"
        argstr = f"{progNameShrt}{mode}_{L}_{p:.3g}_{geo[:3]}_{c[3:]}"
        return slanzarvstr + argstr

else:

    def slanzarv_str(*_):
        return ""


#
if execBool or printBool:
    if execBool and printBool:

        def operate(s, count):
            print(s)
            os.system(s)
            count += 1

    elif execBool:

        def operate(s, count):
            os.system(s)
            count += 1

    elif printBool:

        def operate(s, *args):
            print(s)

    def exec_str(L, p, geo, cell, navg, mode, eigmode, T, bins_count, howmany):
        lnchStr = f"python src/{progName}.py"
        argstr = f"""{L} {p:.3g} -g {geo} -c {cell} -n {navg} --mode={mode} --eigen_mode={eigmode} --period={T} --bins_count={bins_count} --howmany={howmany}"""
        return f"{slanzarv_str(mode, L, p, geo, cell)} {lnchStr} {argstr}"

else:
    exit(0)

count = 0
for L in List:
    for p in plist:
        estring = exec_str(L, p, geo, cell, navg, progMode, eigmode, T, 
                            bins_count, howmany)
        operate(estring, count)