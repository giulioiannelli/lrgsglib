from lrgsglib.config.progargs import *
from parsers.shared import GraphArgumentParser, parse_arguments

# Build combined argument dictionary
# Note: Exclude verbose from action_args_dict to avoid conflict with -v for variant
# We'll add verbose manually with only long form below
common_action_args = {k: v for k, v in action_args_dict.items()
                      if k != tuple(['-v', '--verbose'])}

optionalaction_args_dict = {
    **MC_opt_args,  # MultiplicativeCascade optional args
    **MCG_SlaplSpect_optional_args_dict,  # Spectral analysis optional args
    **MCG_SlaplSpect_action_args_dict,  # Action flags
    **common_action_args,  # Common action flags (print_chrono, etc., excluding verbose)
    # Add verbose with only long form to avoid -v conflict with variant
    tuple(['--verbose']): {
        'help': phelp_verbose,
        'action': argparse.BooleanOptionalAction,
        'default': DEFAULT_VERBOSE
    }
}

# Create parser
parser = GraphArgumentParser(
    description=MCG_SlaplSpect_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Add required positional arguments
for k, v in MCG_SlaplSpect_args.items():
    parser.add_argument(k, **v)

# Add optional/flag arguments
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
# graph-tool eigenvectors are stored node-major (N, k); eigvec_dist reads rows
parser.refuse_engine("gt", "-m eigvec_dist reads NetworkX-layout eigenvectors",
                     when=lambda a: a.mode == "eigvec_dist")
