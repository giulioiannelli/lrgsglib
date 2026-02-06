from .shared import *
from .core import *
from .plotlib import *

sys.setrecursionlimit(DEFAULT_RECURSION_LIMIT)
warnings.simplefilter(action="ignore", category=FutureWarning)

# Note: Logging is opt-in. Use lrgsglib.loglib.setup_custom_logger() if needed.