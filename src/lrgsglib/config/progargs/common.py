from parsers.shared import *
## program helpers
from .phelp.generic import *
## default values
from .defs.generic import *
## argparse tools
from .tools import *
# Program arguments
action_args_dict = {
    tuple(['-rf', '--remove_files']): {
        'help': phelp_remove_files,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_REMOVE_FILES},
    tuple(['-v', '--verbose']): {
        'help': phelp_verbose,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VERBOSE},
    tuple(['-pc', '--print_chrono']): {
        'help': phelp_print_chrono,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_PRINT_CHRONO}
}
# serializer program arguments
srun_opt_args = {
    tuple(['-mMB', '--slanzarv_minMB']): {
        'help': phelp_mMB,
        'type': int,
        'default': DEFAULT_mMB
    },
    tuple(['-MMB', '--slanzarv_maxMB']): {
        'help': phelp_MMB,
        'type': int,
        'default': DEFAULT_MMB
    },
    tuple(['--moretime']): {
        'help': phelp_moretime,
        'type': int,
        'default': DEFAULT_MORETIME
    },
    tuple(['--slanzarv_id']): {
        'help': phelp_slanzarv_id,
        'type': str,
        'default': DEFAULT_SLANZARV_ID
    },
}
# Serializers program arguments
srun_action_args = {
    tuple(['--exec']): {
        'help': phelp_exc,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_EXEC
    },
    tuple(['--print']): {
        'help': phelp_print,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_PRINT
    },
    tuple(['--nomail']): {
        'help': phelp_nomail,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_NOMAIL
    },
    tuple(['--short']): {
        'help': phelp_short,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_SHORT
    },
}