from parsers.ER_IsingDynamics_Serialiser import *
#
args = parser.parse_args()
#
K = args.average_degree
cell = args.cell_type
navg = args.number_of_averages
ic = args.init_cond
outSuffix = args.out_suffix
inSuffix = args.in_suffix
NoClust = args.NoClust
runlang = args.runlang
thrmsteps = args.thermsteps
ic_gs = args.init_cond.startswith('ground_state')
number = int(args.init_cond.split('_')[-1]) if ic_gs else 0
#
progName = ER_IsingDynamics_progName
progNameShrt = ER_IsingDynamics_progNameShrt
execBool = args.exec
printBool = args.print
#
N_list = [1024, 4096]
pflip_list = np.linspace(0.01, 0.8, num=25)
Tlist = np.concatenate([
    np.linspace(0.4, 10, 40),
])

if args.slanzarv_minMB == args.slanzarv_maxMB:
    def memoryfunc(*_):
        return args.slanzarv_minMB
else:
    def memoryfunc(x):
        return int(
            np.interp(
                x,
                [min(N_list), max(N_list)],
                [args.slanzarv_minMB, args.slanzarv_maxMB],
            )
        )
def slanzarv_str(mode, L, p, pflip, T, c):
    slanzarvopt = "--nomail --jobname "
    slanzarvstr = f"slanzarv -m {memoryfunc(L)} {slanzarvopt}"
    argstr = '_'.join([progNameShrt, mode, f"{L}", f"{p:.3g}", f"{pflip:.3g}", 
                       f"{T:.3g}", c, str(number)])
    return slanzarvstr + argstr
#
if execBool or printBool:
    count = 0
    count_exe = 0
    def operate(s):
        global count, count_exe
        if execBool and printBool:
            print(s)
            os.system(s)
            count += 1
            count_exe += 1
        elif execBool:
            os.system(s)
            count_exe += 1
        elif printBool:
            print(s)
            count += 1
    def exec_str(N, p, pflip, T, cell, navg, runlang, inSuffix, out_suffix, 
                 NoClust):
        lnchStr = f"python src/{progName}.py"
        inSuffix = inSuffix or f"p={pflip:.3g}"
        # outStr = f"-os {out_suffix} " if out_suffix else " "
        argStr = " ".join(filter(None, [
            f"{N} {p:.3g} {pflip:.3g} {T:.3g}",
            f"-c {cell} -n {navg} -ic {ic}",
            f"-rl {runlang} -is {inSuffix}",
            f"-os {out_suffix}" if out_suffix else "",
            f"-nc {NoClust}",
            f"-ts {thrmsteps}"
        ])).strip()
        # (f"{L} {p:.3g} {pflip:.3g} {T:.3g} -c {cell} -n {navg} -ic {ic} "
        #           f"-rl {runlang} -is {inSuffix} "
        #           f"{outStr}"
        #           f"-nc {NoClust}")
        slanzStr = slanzarv_str(runlang, N, p, pflip, T, cell)
        return ' '.join([slanzStr, lnchStr, argStr])

else:
    exit(0)

for N in N_list:
    p = K/N
    for pflip in pflip_list:
            for T in Tlist:
                estring = exec_str(N, p, pflip, T, cell, navg, 
                                   runlang, inSuffix, outSuffix, NoClust)
                operate(estring)
print(f"Total number of jobs executed: {count_exe}")
print(f"Total number of jobs printed: {count}")