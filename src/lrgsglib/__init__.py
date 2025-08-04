from .shared import *
from .core import *
from .plotlib import *
from .loglib import setup_custom_logger

sys.setrecursionlimit(DEFAULT_RECURSION_LIMIT)
warnings.simplefilter(action="ignore", category=FutureWarning)

# Automatically configure on import (you can also require an explicit call)
setup_custom_logger(__name__)