from __future__ import annotations

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
import numpy as np


ext_modules = [
    Pybind11Extension(
        "lrgsglib._random_walk",
        ["src/lrgsglib/Ccore/statsys/random_walk.cpp"],
        include_dirs=[np.get_include()],
        cxx_std=17,
        extra_compile_args=["-O3"],
    ),
]

setup(cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)
