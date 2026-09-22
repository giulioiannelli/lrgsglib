#
import argparse
import copy
import dis
import glob
import itertools
import operator
import os
import random
import re
import string
import struct
import subprocess
import sys
import time
import warnings

import lmfit
import powerlaw
import scipy

#
try:  # cupy is GPU-optional: importing lrgsglib must not require CUDA.
    import cupy as cp
except Exception:  # ImportError, or CUDA/driver errors raised at import time
    cp = None
import pickle as pk

#
from collections import Counter
from collections.abc import Iterable
from decimal import Decimal
from fractions import Fraction
from itertools import product
from numbers import Number
from operator import itemgetter
from os.path import join as pth_join
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    Type,
    Union,
)

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from cycler import cycler
from joblib import Memory
from networkx.classes.graph import Graph
from networkx.drawing.layout import rescale_layout
from numpy.typing import NDArray
from PIL import Image
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import (
    cophenet,
    dendrogram,
    fcluster,
    leaves_list,
    linkage,
)
from scipy.interpolate import griddata, pchip
from scipy.io import loadmat
from scipy.linalg import eigvalsh as seigvalsh
from scipy.linalg import expm, fractional_matrix_power
from scipy.ndimage import gaussian_filter, gaussian_filter1d, zoom
from scipy.optimize import curve_fit
from scipy.signal import (
    argrelextrema,
    butter,
    find_peaks,
    medfilt,
    peak_prominences,
    sosfiltfilt,
)
from scipy.sparse import coo_matrix, csr_array, csr_matrix, diags
from scipy.sparse import identity as scsp_identity
from scipy.sparse import spdiags
from scipy.sparse.linalg import eigsh as scsp_eigsh
from scipy.sparse.linalg import expm as sparse_expm
from scipy.spatial.distance import pdist, squareform
from scipy.stats import gaussian_kde
from sklearn.datasets import fetch_openml
from tqdm.auto import tqdm  # auto-detects Jupyter vs terminal for clean output

# No `__all__`: this module is the canonical kitchen-sink import for
# notebooks/labs (`from lrgsglib.notebooks import *` chains through here).
# Top-level package surface is curated explicitly in `lrgsglib/__init__.py`.
