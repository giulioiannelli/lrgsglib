"""Parser for L3D_TopoComparison_Serializer.

Sweeps over side lengths and pflip values. Temperature is only relevant
for SA/topo_fca (annealing schedules), not for topo_met.  The comparison
serialiser fixes the three runlangs internally.
"""

from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments

TOPO_COMP_description = """Batch comparison serialiser: SA vs topo_met vs topo_fca on 3D lattices."""

# Sweep over (side, pflip) only — temperature handled by SA schedule
_sweep_args = {
    tuple(['-s1', '--side1_list']): {
        'help': phelp_side1_list,
        'type': int,
        'nargs': '+',
        'default': DEFAULT_SIDE1_LIST,
    },
    tuple(['-pFT', '--pflip_linsp']): {
        'help': phelp_pflip_linsp,
        'type': parse_multiple_linspace,
        'default': DEFAULT_PFLIP_LINSP,
    },
}

optionalaction_args_dict = {
    **_sweep_args,
    **L3D_opt_args,                 # geometry, pdil, etc.
    **IsDyn_sa_args,                # SA schedule (shared by SA and topo_fca)
    **IsDyn_topo_args,              # topo params
    **IsDyn_srun_args,              # slanzarv memory, id, etc.
    **L3D_ISDYN_srun_action_args,   # exec, print, nomail, short
    # Additional overrides
    tuple(['-na', '--number_of_averages']): {
        'help': 'Number of disorder realizations per method',
        'type': int,
        'default': 100,
    },
    tuple(['--n_thermal']): {
        'help': phelp_n_thermal,
        'type': int,
        'default': DEFAULT_N_THERMAL,
    },
    tuple(['-ic', '--init_cond']): {
        'help': phelp_ic,
        'type': str,
        'default': 'rand',
    },
    tuple(['-fq', '--freq']): {
        'help': phelp_freq,
        'type': int,
        'default': DEFAULT_FREQ,
    },
    tuple(['-ts', '--thrmsteps']): {
        'help': phelp_thrmsteps,
        'type': int,
        'default': DEFAULT_THRMSTEPS,
    },
}

parser = argparse.ArgumentParser(
    description=TOPO_COMP_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
)

for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
