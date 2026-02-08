from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import distutils.sysconfig as sysconfig
import os
import sys

class get_boost_include(object):
    """Helper class to determine the Boost include path"""
    def __str__(self):
        # Path to Boost headers in the Conda environment
        return os.path.join(env_path, 'include')

def get_python_lib():
    """Helper function to determine the Python library path"""
    return sysconfig.get_python_lib()

def get_python_version():
    """Helper function to get the Python version in the format Boost expects"""
    return f"{sys.version_info.major}{sys.version_info.minor}"

# Get the environment path
env_path = os.path.dirname(sys.executable)

ext_modules = [
    Extension(
        'ising_model_store',
        ['main.cpp', 'ising_model_store.cpp'],  # Include all relevant source files
        include_dirs=[
            # Path to Boost headers in the Conda environment
            get_boost_include(),
            # Include dirs for the Conda environment
            os.path.join(env_path, 'include'),
            # Additional include directories if needed
            # os.path.join(env_path, 'include', 'boost'),
        ],
        library_dirs=[
            # Path to Python library
            get_python_lib(),
            # Library dirs for the Conda environment
            os.path.join(env_path, 'lib'),
        ],
        extra_compile_args=['-O3'],  # Enable O3 optimization
        extra_link_args=['-O3'],     # Enable O3 optimization for linking
        libraries=[f'boost_python{get_python_version()}'],  # Dynamically set the Boost.Python library
        language='c++'
    ),
]

setup(
    name='ising_model_store',
    version='0.0.1',
    author='Your Name',
    author_email='your.email@example.com',
    description='An Ising model simulation using Boost.Python',
    ext_modules=ext_modules,
    install_requires=['boost'],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
