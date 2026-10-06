DEFAULT_BINSC = 500
DEFAULT_EIGMODE = "smallest"  # Eigenvector selection mode
DEFAULT_EXEC = False
DEFAULT_FREQ = 2
DEFAULT_INSFFX = ""
DEFAULT_mMB = 2**11
DEFAULT_MMB = 2**11
DEFAULT_MORETIME = 0
DEFAULT_NAVG = 500
DEFAULT_NOMAIL = True
DEFAULT_OUTSFFX = ""
DEFAULT_PRINT = False
DEFAULT_PRINT_CHRONO = False
DEFAULT_REMOVE_FILES = True
DEFAULT_SAVING_FREQUENCY = 10
DEFAULT_SHORT = False
DEFAULT_SLANZARV_ID = ""
DEFAULT_GPU = False
DEFAULT_GPUTYPE = None
DEFAULT_PARTITION = None
DEFAULT_VAL = "=1"
DEFAULT_VERBOSE = False
DEFAULT_WORKDIR = ""
DEFAULT_GRAPH_ENGINE = "nx"
GRAPH_ENGINE_CHOICES = ("nx", "gt")
## default values depending on default arguments
DEFAULT_DATA_SAVE_FREQ = DEFAULT_NAVG // 10
# BLAS/OpenMP threads per array task (sbatch exports the same caps)
DEFAULT_SRUN_THREADS = 1
# SLURM arrays (kernels/Serializer.render_sbatch_tiles_script)
DEFAULT_ARRAY_MAX_SIZE = 1000  # tasks per array; keep below MaxArraySize
SRUN_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
SBATCH_ENV = {"PYTHONUNBUFFERED": "1"}  # exported in every array script
SBATCH_LOG_DIR = ".log"  # task logs, relative to the job's --chdir
# #SBATCH values are not shell-quoted: sbatch splits them at whitespace, ends
# them at "#", and strips quotes and backslashes, so these are refused.
SBATCH_FORBIDDEN_CHARS = "#'\"\\"
SBATCH_ENV_PREFIX = "SBATCH_"  # input variables that override #SBATCH lines
