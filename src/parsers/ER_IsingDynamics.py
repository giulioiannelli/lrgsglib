from lrgsglib.core import *
from parsers.shared import *
#
def parse_fraction(value):
    try:
        if '/' in value:
            return float(Fraction(value))
        else:
            return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid number: {value}")
#
description = """
    Computational resourses regarding the Ising Dynamics on signed Erdos-Renyi 
    graphs.
"""
phelp_N = """
    Number of nodes of the graph.
"""
phelp_p = """
    Erdos-Renyi edge probability.
"""
phelp_pflip = """
    Edge flipping probability.
"""
phelp_T = """
    Temperature of the Ising model.
"""
#
DEFAULT_INIT_COND = 'ground_state_0'
DEFAULT_CELL = 'rand'
DEFAULT_RUNLANG = 'C1'
DEFAULT_NAVG = 500
DEFAULT_NAVG2 = 10
DEFAULT_INSFFX = ''
DEFAULT_OUTSFFX = ''
DEFAULT_NOCLUST = 1
DEFAULT_REMOVE_FILES = True
DEFAULT_WORKDIR = ''
DEFAULT_MAX_ITER_ER_GC = 20
DEFAULT_PRINT_CHRONO = False
DEFAULT_THRMSTEPS = 20
DEFAULT_VERBOSE = False
#
phelp_thrmsteps = """
    Number of thermalization steps.
"""
phelp_remove_files = """
    Remove the input files after the computation.
"""
phelp_print_chrono = """
    Print the chronometer.
"""
phelp_cell = f"""
    Topological defect class: 'rand', 'randXERR'.
"""
phelp_navg = f"""
    Number of averages to compute.
"""
phelp_initCond = f"""
    Initial condition for the Ising model.
"""
phelp_runlang = f"""
    Language for running the Ising model
"""
phelp_inSuffix = f"""
    Suffix for the input files
"""
phelp_outSuffix = f"""
    Suffix for the output files
"""
phelp_NoClust = f"""
    Number of clusters to compute
"""
phelp_workdir = f"""
    Working directory
"""
phelp_navg2 = f"""
    Number of averages for the cluster statistic
"""
phelp_verbose = f"""
    Verbose mode
"""
#
parsDict = {
    'N': {'help': phelp_N, 'type': int},
    'p': {'help': phelp_p, 'type': parse_fraction},
    'pflip' : {'help': phelp_pflip, 'type': float},
    'T': {'help': phelp_T, 'type': float}
}
#
parsDictOpt = {
    'cell': {'names': ['-c', '--cell_type'],
            'help': phelp_cell,
            'type': str,
            'default': DEFAULT_CELL},
    'navg': {'names': ['-n', '--number_of_averages'],
            'help': phelp_navg,
            'type': int,
            'default': DEFAULT_NAVG},
    'init_cond': {'names': ['-ic', '--init_cond'],
                    'help': phelp_initCond,
                    'type': str,
                    'default': DEFAULT_INIT_COND},
    'runlang': {'names': ['-rl', '--runlang'],
                'help': phelp_runlang,
                'type': str,
                'default': DEFAULT_RUNLANG},
    'out_suffix': {'names': ['-os', '--out_suffix'],
                    'help': phelp_outSuffix,
                    'type': str,
                    'default': DEFAULT_OUTSFFX},
    'NoClust': {'names': ['-nc', '--NoClust'],
                'help': phelp_NoClust,
                'type': int,
                'default': DEFAULT_NOCLUST},
    'workdir': {'names': ['-wd', '--workdir'],
                'help': phelp_workdir,
                'type': str,
                'default': DEFAULT_WORKDIR},
    'navg2': {'names': ['-n2', '--number_of_averages2'],
                'help': phelp_navg2,
                'type': int,
                'default': DEFAULT_NAVG2},
    'thrmsteps': {'names': ['-ts', '--thrmsteps'],
                'help': phelp_thrmsteps,
                'type': int,
                'default': DEFAULT_THRMSTEPS}
}
#
parDA = {'remove_files': {'names': ['-rf', '--remove_files'],
                            'help': phelp_remove_files,
                            'action': argparse.BooleanOptionalAction,
                            'default': DEFAULT_REMOVE_FILES},
        'print_chrono': {'names': ['-pc', '--print_chrono'],
                        'help': phelp_print_chrono,
                        'action': argparse.BooleanOptionalAction,
                        'default': DEFAULT_PRINT_CHRONO},
        'verbose': {'names': ['-v', '--verbose'],
                    'help': phelp_verbose,
                    'action': argparse.BooleanOptionalAction,
                    'default': DEFAULT_VERBOSE}
}
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# Mandatory arguments
def parsDict_get(var, key):
    return parsDict[var].get(key, None)
for v in parsDict.keys():
    parser.add_argument(v,
        help=parsDict_get(v, 'help'),
        type=parsDict_get(v, 'type'),
    )
# Optional parameters
def parsDict_get(var, key):
    return parsDictOpt[var].get(key, None)
for ov in parsDictOpt.keys():
    parser.add_argument(*parsDict_get(ov, 'names'),
        default=parsDict_get(ov, 'default'),
        help=parsDict_get(ov, 'help'),
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