from parsers.L3D_TransCluster_Serialiser import *
#
args = parser.parse_args()
#
progName = L3D_TransCluster_progName
progNameShrt = L3D_TransCluster_progNameShrt
navg = args.number_of_averages
fullMode = args.mode
progMode = args.mode.split('_')[-1]
execBool = args.exec
printBool = args.print
typf = args.float_type
outsx = args.out_suffix
pdil  = args.pdil
mu = args.mu
edge_weight = args.edge_weight
#
match fullMode:
    case s if s.endswith('pCluster'):
        List = [10, 20, 30]
        geo = args.geometry
        cell = args.cell_type
        if mu != parser.get_default('mu'):
            siglist = np.linspace(0, 1.5 * mu, 100)
        else:
            plist = np.linspace(0.06, .3, num=20)
    case s if s.endswith('ordParam'):
        List = [10, 20, 30, 40]
        def linspacepfunc(*_):
            return np.linspace(0, 0.5, 100)
        geometry_cell_dict = {'simple_cubic': ['rand']}
        if mu != parser.get_default('mu'):
            siglist = {geo: {cell: np.linspace(0, 1.5 * mu, 100) for cell in cells}  
                    for geo, cells in geometry_cell_dict.items()}
        else:
            plist = {geo: {cell: linspacepfunc(geo, cell) for cell in cells}  
                    for geo, cells in geometry_cell_dict.items()}
    case _:
        raise ValueError("Invalid mode specified")
#
match fullMode:
    case s if s.startswith('slanzarv'):
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
    case _:
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
    def exec_str(L, p, geo, cell, sigma=0, navg=navg, mode=progMode, pdil=pdil):
        lnchStr = f"python src/{progName}.py"
        optStr = f"-o {outsx}" if outsx else ""
        argstr = (f"{L} {p:.3g} -g {geo} -c {cell} -n {navg} -t {typf} {optStr}"+
                  f" --mode={mode}" + f" --pdil={pdil:.3g}" + f" --mu={mu:.3g}" + 
                  f" --sigma={sigma:.3g}" + f" --edge_weight={edge_weight}")
        return (f"{slanzarv_str(mode, L, p, geo, cell)} "
                        f"{lnchStr} {argstr}")
else:
    exit(0)

count = 0
if fullMode.endswith('pCluster'):
    for L in List:
        match edge_weight:
            case 'normal':
                for sigma in siglist:
                    estring = exec_str(L, 0, geo, cell, sigma=sigma)
                    count = operate(estring, count)
            case 'flip':
                for p in plist:
                    estring = exec_str(L, p, geo, cell)
                    count = operate(estring, count)
elif fullMode.endswith('ordParam'):
    for L in List:
        for geo, cellst in geometry_cell_dict.items():
            for cell in cellst:
                match edge_weight:
                    case 'normal':
                        for sigma in siglist[geo][cell]:
                            estring = exec_str(L, 0, geo, cell, sigma=sigma)
                            count = operate(estring, count)
                    case 'flip':
                        for p in plist[geo][cell]:
                            estring = exec_str(L, p, geo, cell)
                            count = operate(estring, count)
else:
    print(f"no program executed, unknown mode provided: {fullMode}")
    exit(0)         
print("submitted jobs: ", count)