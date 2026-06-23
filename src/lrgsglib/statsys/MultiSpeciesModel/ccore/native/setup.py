"""Build script for the _multispec_native pybind11 extension."""

from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext as _build_ext
from setuptools import setup


class build_ext(_build_ext):
    """Strip C++-only flags when compiling plain C sources."""

    def build_extensions(self):
        original_compile = self.compiler._compile

        def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
            if src.endswith(".c"):
                extra_postargs = [
                    f for f in extra_postargs
                    if not f.startswith(("-std=c++", "-std=gnu++"))
                ]
            original_compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

        self.compiler._compile = _compile
        super().build_extensions()


HERE = Path(__file__).resolve().parent
CCORE_DIR = HERE.parent                              # MultiSpeciesModel/ccore
SHARED_CCORE = CCORE_DIR.parent.parent / "_ccore"    # statsys/_ccore
SFMT_DIR = SHARED_CCORE / "SFMT"
# GCC 15 + old conda sysroot compat (see build/gcc15_compat.h)
LRGSG_ROOT = CCORE_DIR.parents[4]  # ccore -> MultiSpeciesModel -> statsys -> lrgsglib -> src -> ROOT
GCC15_COMPAT = str(LRGSG_ROOT / "build" / "gcc15_compat.h")

ext_modules = [
    Pybind11Extension(
        "_multispec_native",
        sources=[
            str(CCORE_DIR / "multispec_native.cpp"),
            str(CCORE_DIR / "LRGSG_multispec.c"),
            str(SHARED_CCORE / "LRGSG_vecdynsys.c"),
            str(SHARED_CCORE / "LRGSG_utils.c"),
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
            ("SFMT_MEXP", "19937"),
            ("HAVE_SSE2", None),
        ],
        extra_compile_args=["-O3", "-msse2", "-include", GCC15_COMPAT],
    ),
]

setup(
    name="_multispec_native",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
