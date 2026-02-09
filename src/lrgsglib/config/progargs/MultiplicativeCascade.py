from .common import *

# Import definitions and help strings
from .defs.MultiplicativeCascade import *
from .phelp.MultiplicativeCascade import *

# Program metadata
MC_progName = 'MultiplicativeCascade'
MCG_MCGSS_progname = 'MCG_SlaplSpect'
MCG_MCGSS_progname_shrt = 'MCGSS'

# Argument dictionary for MC graphs (positional arguments)
MC_args = {
    'p1': {
        'help': phelp_mc_p1,
        'type': float,
    },
    'p2': {
        'help': phelp_mc_p2,
        'type': float,
    },
    'p3': {
        'help': phelp_mc_p3,
        'type': float,
    },
    'p4': {
        'help': phelp_mc_p4,
        'type': float,
    },
    'fraction': {
        'help': phelp_mc_fraction,
        'type': float,
    },
    'iterations': {
        'help': phelp_mc_iterations,
        'type': int,
    }
}

# Optional arguments (flags/switches)
MC_opt_args = {
    tuple(['-st', '--stochastic']): {
        'help': phelp_mc_stochastic,
        'action': 'store_true',
        'default': DEFAULT_MC_STOCHASTIC
    },
    tuple(['-pbc', '--periodic']): {
        'help': phelp_mc_periodic,
        'action': 'store_true',
        'default': DEFAULT_MC_PERIODIC
    },
    tuple(['-v', '--variant']): {
        'help': phelp_mc_variant,
        'type': str,
        'choices': ['standard', 'exp_clocks'],
        'default': DEFAULT_MC_VARIANT
    },
    tuple(['-p', '--pflip']): {
        'help': phelp_p,  # Reuse from SignedGraph (probability of flipping edge signs)
        'type': float,
        'default': 0.0
    },
    tuple(['-na', '--number_of_averages']): {
        'help': phelp_navg,  # Reuse from common
        'type': int,
        'default': DEFAULT_NAVG
    },
    tuple(['-os', '--out_suffix']): {
        'help': 'Optional suffix appended to output directory name',
        'type': str,
        'default': ""
    },
    tuple(['-wd', '--workdir']): {
        'help': phelp_workdir,  # Reuse from common
        'type': str,
        'default': DEFAULT_WORKDIR
    },
    tuple(['-ge', '--graph_engine']): {
        'help': phelp_graph_engine,
        'type': str,
        'default': DEFAULT_GRAPH_ENGINE,
        'choices': ['nx', 'gt'],
    },
}

# Serializer-specific argument lists (for parameter sweeps)
MCG_SlaplSpect_srun_description = f"""
    Batch launcher for {MCG_MCGSS_progname}.py:
    Compute signed Laplacian spectra for MultiplicativeCascadeGraph over parameter sweeps.
"""

MCG_SlaplSpect_srun_optional_args_dict = {
    # Cascade probability matrix specification
    tuple(['--matrix_list']): {
        'help': 'List of 2x2 cascade matrices as semicolon-separated 4-tuples: "p1,p2,p3,p4;p1,p2,p3,p4;..."',
        'type': str,
        'default': None
    },
    tuple(['--matrix_file']): {
        'help': 'File containing matrices (one per line: p1 p2 p3 p4)',
        'type': str,
        'default': None
    },
    # fraction sweep
    tuple(['--fraction_list']): {
        'help': 'List of fraction values (comma-separated)',
        'type': lambda s: [float(x.strip()) for x in s.split(',')],
        'default': None
    },
    tuple(['--fraction_linsp']): {
        'help': 'Linear space for fraction: start,stop,num',
        'type': lambda s: [float(x.strip()) for x in s.split(',')],
        'default': None
    },
    # iterations sweep
    tuple(['--iterations_list']): {
        'help': 'List of iteration values (comma-separated)',
        'type': lambda s: [int(x.strip()) for x in s.split(',')],
        'default': None
    },
    tuple(['--iterations_linsp']): {
        'help': 'Linear space for iterations: start,stop,num',
        'type': lambda s: [int(x.strip()) for x in s.split(',')],
        'default': None
    },
    # pflip sweep (edge sign flipping probability)
    tuple(['--p_list']): {
        'help': 'List of pflip values (edge sign flip probability)',
        'type': lambda s: [float(x.strip()) for x in s.split(',')],
        'default': None
    },
    tuple(['--p_linsp']): {
        'help': 'Linear space for pflip: start,stop,num',
        'type': lambda s: [float(x.strip()) for x in s.split(',')],
        'default': None
    }
}
