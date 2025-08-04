from lrgsglib.proglib import *
#
parser = argparse.ArgumentParser(
    description=L3D_Recon_srun_description, 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False
)
#
for k,v in {**L3D_Recon_srun_optional_args_dict, **srun_action_args}.items():
    parser.add_argument(*k, **v)