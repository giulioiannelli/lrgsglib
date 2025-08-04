from parsers.L2D_TransCluster_Serialiser import *
#
args = parser.parse_args()
#
progName = Lattice2D_TransCluster_progName
progNameShrt = Lattice2D_TransCluster_progNameShrt
navg = args.number_of_averages
fullMode = args.mode
progMode = args.mode.split('_')[-1]
execBool = args.exec
printBool = args.print
typf = args.float_type
outsx = args.out_suffix
prew = args.prew
#
if fullMode.endswith('pCluster'):
    List = [16, 32, 64, 96, 128]
    geo = args.geometry
    cell = args.cell_type
    plist = np.linspace(0.06, .3, num=25)
elif fullMode.endswith('ordParam'):
    List = [32, 64, 96, 128]
    def linspacepfunc(geo, cell):
        if cell != 'randXERR':
            if geo != 'hexagonal':
                return np.linspace(0, 0.3, 100)
            else:
                return np.linspace(0, 1, 100)
        return np.linspace(0, 1, 100)
    geometry_cell_dict = {'hexagonal': L2D_RAND_CELL_LIST}
    # geometry_cell_dict = {'squared': ['rand'],
    #                       'triangular': ['rand']}
    plist = {geo: {cell: linspacepfunc(geo, cell) for cell in cells}  
            for geo,cells in geometry_cell_dict.items()}
#
if fullMode.startswith('slanzarv'):
    if args.slanzarv_minMB == args.slanzarv_maxMB:
        def memoryfunc(*_):
            return args.slanzarv_minMB
    else:
        def memoryfunc(x):
            return int(np.interp(x, [min(List), max(List)], 
                             [args.slanzarv_minMB, args.slanzarv_maxMB]))
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
            return count +1
    elif execBool:
        def operate(s, count):
            os.system(s)
            return count +1
    elif printBool:
        def operate(s, *args):
            print(s)
            return count +1
    def exec_str(L, p, geo, cell, navg=navg, mode=progMode):
        lnchStr = f"python src/{progName}.py"
        optStr = f"-o {outsx}" if outsx else ""
        argstr = (f"{L} {p:.3g} -g {geo} -c {cell} -n {navg} -t {typf} {optStr} "+
                  f"--mode={mode} --prew={prew}")
        return (f"{slanzarv_str(mode, L, p, geo, cell)} "
                        f"{lnchStr} {argstr}")
else:
    exit(0)

count = 0
if fullMode.endswith('pCluster'):
    for L in List:
        for p in plist:
            estring = exec_str(L, p)
            count = operate(estring, count)
elif fullMode.endswith('ordParam'):
    for L in List:
        for geo, cellst in geometry_cell_dict.items():
            for cell in cellst:
                for p in plist[geo][cell]:
                    estring = exec_str(L, p, geo, cell)
                    count = operate(estring, count)
else:
    print(f"no program executed, unkonw mode provided: {fullMode}")
    exit(0)         
print("submitted jobs: ", count)