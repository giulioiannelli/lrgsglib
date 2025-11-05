from lrgsglib.config.progargs import *

optionalaction_args_dict = {
    **SCSGeneralized_optional_args_dict,
    **SCS_TransCluster_optional_args_dict,
    **SCS_TransCluster_action_args_dict,
}

parser = argparse.ArgumentParser(
    description=SCS_TransCluster_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for k, v in SCS_TransCluster_args.items():
    parser.add_argument(k, **v)
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
