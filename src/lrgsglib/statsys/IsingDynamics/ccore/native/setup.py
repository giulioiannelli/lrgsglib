"""Build script for _ising_native pybind11 extension."""

from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

HERE = Path(__file__).resolve().parent
CCORE_DIR = HERE.parent          # IsingDynamics/ccore
SHARED_CCORE = CCORE_DIR.parent.parent / "_ccore"  # statsys/_ccore
SFMT_DIR = SHARED_CCORE / "SFMT"
# GCC 15 + old conda sysroot compat (see build/gcc15_compat.h)
LRGSG_ROOT = CCORE_DIR.parents[4]  # ccore → IsingDynamics → statsys → lrgsglib(pkg) → src → ROOT
GCC15_COMPAT = str(LRGSG_ROOT / "build" / "gcc15_compat.h")

ext_modules = [
    Pybind11Extension(
        "_ising_native",
        sources=[
            str(CCORE_DIR / "ising_native.cpp"),
            str(CCORE_DIR / "LRGSG_rbim.c"),
            str(CCORE_DIR / "LRGSG_sa.c"),
            str(CCORE_DIR / "LRGSG_pt.c"),
            str(CCORE_DIR / "LRGSG_cluster.c"),
            str(SHARED_CCORE / "LRGSG_utils.c"),
            str(SHARED_CCORE / "LRGSG_uf.c"),
            str(SHARED_CCORE / "sfmtrng.c"),
            str(SFMT_DIR / "SFMT.c"),
            str(SHARED_CCORE / "LRGSG_bindynsys.c"),
        ],
        include_dirs=[
            str(SHARED_CCORE),
            str(SFMT_DIR),
            str(CCORE_DIR),
        ],
        define_macros=[
            ("DSFMT_MEXP", "19937"),
            ("HAVE_SSE2", None),
        ],
        extra_compile_args=["-O3", "-msse2", "-include", GCC15_COMPAT],
    ),
]

setup(
    name="_ising_native",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
