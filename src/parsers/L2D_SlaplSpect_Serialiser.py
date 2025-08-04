from lrgsglib.core import *
#
description = """
    Serialiser for Lattice2D_SlaplSpect.py
"""
L2D_SlaplSpect_progName = "L2D_SlaplSpect"
L2D_SlaplSpect_progNameShrt = "L2DSS"
# Default values for the optional parameters
DEFAULT_BINSC = 500
DEFAULT_EIGMODE = 'scipy'
DEFAULT_NAVG = 10**4
DEFAULT_PERIOD = DEFAULT_NAVG//20
DEFAULT_WORKDIR = ''
DEFAULT_HOWMANY = 1
DEFAULT_GEO = 'squared'
DEFAULT_CELL = 'rand'
DEFAULT_MODE = 'slanzarv_eigvec_dist'
#
DEFAULT_PRINT = False
DEFAULT_EXEC = False
DEFAULT_mMB = 2**10
DEFAULT_MMB = 2**14
# Helpers for argparse descriptions
phelp_print = f"""
    Option to print the output of the Serialiser. 
"""
phelp_exc = f"""
    Option to exec the output of the Serialiser.
"""
phelp_mMB = f"""
    Minimum MB quantity to be allocated for the single process 
"""
phelp_MMB = f"""
    Maximum MB quantity to be allocated for the single process 
"""
phelp_binsc = f"""
    Number of bins for the distribution sampling 
"""
phelp_cell = f"""
    Topological defect class: 'rand', 'randXERR', 'randZERR', 'ball_<R>' 
    with type(<R>)=int. 
"""
phelp_eigMode = f"""
    Spectral computations (numpy/scipy) 
"""
phelp_geo = f"""
    Geometry of the lattice. 
"""
phelp_howmany = f"""
    Number of eigenvalues to compute 
"""
phelp_mode = f"""
    Mode, (eigDistr, slanzarv_eigDistr) 
"""
phelp_navg = f"""
    Number of averages to compute 
"""
phelp_period = f"""
    Period for saving the data 
"""
phelp_workDir = f"""
    Working directory 
"""
#
parDO = {'mode': {'names': ['-m', '--mode'], 
                'help': phelp_mode, 
                'type': str, 
                'default': DEFAULT_MODE},
        'eigmode': {'names': ['-em', '--eigen_mode'], 
                'help': phelp_eigMode, 
                'type': str, 
                'default': DEFAULT_EIGMODE},
        'geo': {'names': ['-g', '--geometry'],
                'help': phelp_geo,
                'type': str,
                'default': DEFAULT_GEO} ,
        'cell': {'names': ['-c', '--cell_type'],
                'help': phelp_cell,
                'type': str,
                'default': DEFAULT_CELL},
        'navg': {'names': ['-n', '--number_of_averages'],
                'help': phelp_navg,
                'type': int,
                'default': DEFAULT_NAVG},
        'period': {'names': ['-prd', '--period'],
                'help': phelp_period,
                'type': int,
                'default': DEFAULT_PERIOD},
        'bins_count': {'names': ['-bc', '--bins_count'],
                'help': phelp_binsc,
                'type': int,
                'default': DEFAULT_BINSC},
        'howmany': {'names': ['-hm', '--howmany'],
                'help': phelp_howmany,
                'type': int,
                'default': DEFAULT_HOWMANY},
        'workDir': {'names': ['-wd', '--workDir'],
                'help': phelp_workDir,
                'type': str,
                'default': DEFAULT_WORKDIR},
        'mMB': {'names': ['-mMB', '--slanzarv_minMB'],               
                'help': phelp_mMB,
                'type': int,
                'default': DEFAULT_mMB},
        'MMB': {'names': ['-MMB', '--slanzarv_maxMB'],
                'help': phelp_MMB,
                'type': int,
                'default': DEFAULT_MMB}
}
parDA = {'exec': {'names': ['-e', '--exec'],
                        'help': phelp_exc,
                        'action': argparse.BooleanOptionalAction,
                        'default': DEFAULT_EXEC},
                'print': {'names': ['-p', '--print'],                 
                        'help': phelp_print,
                        'action': argparse.BooleanOptionalAction,
                        'default': DEFAULT_PRINT}
                        }

# Setup the argument parser
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Optional parameters
def parsDict_get(var, key):
    return parDO[var].get(key, None)
for ov in parDO.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help')+\
            f"(default: {parsDict_get(ov, 'default')})",
        type=parsDict_get(ov, 'type'),
    )
# Optional parameters
def parsDict_get(var, key):
    return parDA[var].get(key, None)
for ov in parDA.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help'),
        action=parsDict_get(ov, 'action')
    )