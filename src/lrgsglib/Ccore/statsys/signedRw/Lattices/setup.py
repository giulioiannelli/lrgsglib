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
        'random_walker',
        ['random_walker.cpp'],
        include_dirs=[
            # Path to Boost headers in the Conda environment
            get_boost_include(),
            # Include dirs for the Conda environment
            os.path.join(env_path, 'include'),
        ],
        library_dirs=[
            # Path to Python library
            get_python_lib(),
            # Library dirs for the Conda environment
            os.path.join(env_path, 'lib'),
        ],
        libraries=[f'boost_python{get_python_version()}'],  # Dynamically set the Boost.Python library
        language='c++'
    ),
]

setup(
    name='random_walker',
    version='0.0.1',
    author='Your Name',
    author_email='your.email@example.com',
    description='A random walker simulation using Boost.Python',
    ext_modules=ext_modules,
    install_requires=['boost'],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
