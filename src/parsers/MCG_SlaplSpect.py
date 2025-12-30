from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments

# Build combined argument dictionary
optionalaction_args_dict = {
    **MC_opt_args,  # MultiplicativeCascade optional args
    **MCG_SlaplSpect_optional_args_dict,  # Spectral analysis optional args
    **MCG_SlaplSpect_action_args_dict  # Action flags
}

# Create parser
parser = argparse.ArgumentParser(
    description=MCG_SlaplSpect_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Add required positional arguments
for k, v in MCG_SlaplSpect_args.items():
    parser.add_argument(k, **v)

# Add optional/flag arguments
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
