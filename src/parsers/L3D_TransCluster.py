from lrgsglib.config.progargs import *

optionalaction_args_dict = {
    **L3D_TransCluster_optional_args_dict,
    **L3D_TransCluster_action_args_dict,
}

parser = argparse.ArgumentParser(
    description=L3D_TransCluster_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for k, v in L3D_TransCluster_args.items():
    parser.add_argument(k, **v)
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
