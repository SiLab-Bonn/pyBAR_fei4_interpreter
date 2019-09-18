#!/usr/bin/env python

# Building the C++ libraries:
# python setup.py build_ext
#
# Installing package from sources:
# python setup.py install
# For developers (creating a link to the sources):
# python setup.py develop

from setuptools import setup, Extension, find_packages
from distutils.command.build_ext import build_ext
from Cython.Build import cythonize
import numpy as np
import os

copt = {'msvc': ['-Ipybar_fei4_interpreter/external', '/EHsc']}  # Set additional include path and EHsc exception handling for VS
lopt = {}


class build_ext_opt(build_ext):
    def initialize_options(self):
        build_ext.initialize_options(self)
        self.compiler = 'msvc' if os.name == 'nt' else None  # in Anaconda the libpython package includes the MinGW import libraries and a file (Lib/distutils/distutils.cfg) which sets the default compiler to mingw32. Alternatively try conda remove libpython.

    def build_extensions(self):
        c = self.compiler.compiler_type
        if c in copt:
            for e in self.extensions:
                e.extra_compile_args = copt[c]
        if c in lopt:
            for e in self.extensions:
                e.extra_link_args = lopt[c]
        build_ext.build_extensions(self)


extensions = [
    Extension('pybar_fei4_interpreter.data_interpreter', ['pybar_fei4_interpreter/data_interpreter.pyx', 'pybar_fei4_interpreter/Interpret.cpp', 'pybar_fei4_interpreter/Basis.cpp']),
    Extension('pybar_fei4_interpreter.data_histograming', ['pybar_fei4_interpreter/data_histograming.pyx', 'pybar_fei4_interpreter/Histogram.cpp', 'pybar_fei4_interpreter/Basis.cpp']),
    Extension('pybar_fei4_interpreter.analysis_functions', ['pybar_fei4_interpreter/analysis_functions.pyx'])
]


f = open('VERSION', 'r')
version = f.readline().strip()
f.close()

author = 'Jens Janssen, David-Leon Pohl'
author_email = 'janssen@physik.uni-bonn.de, pohl@physik.uni-bonn.de'

# Requirements from requirements.txt
with open('requirements.txt') as f:
    install_requires = f.read().splitlines()

setup(
    name='pyBAR_fei4_interpreter',
    version=version,
    description='pyBAR_fei4_interpreter - Fast ATLAS FE-I4 raw data interpreter for pyBAR',
    url='https://github.com/SiLab-Bonn/pyBAR_fei4_interpreter',
    license='BSD 3-Clause ("BSD New" or "BSD Simplified") License',
    long_description='Interpreter for ATLAS FE-I4A/B raw data for the readout framework pyBAR. It also provides histogramming functions. The interpreter is written in C++ to achieve a high throughput.',
    author=author,
    maintainer=author,
    author_email=author_email,
    maintainer_email=author_email,
    install_requires=install_requires,
    packages=find_packages(),  # exclude=['*.tests', '*.test']),
    include_package_data=True,  # accept all data files and directories matched by MANIFEST.in or found in source control
    package_data={'': ['README.*', 'VERSION'], 'docs': ['*'], 'examples': ['*']},
    ext_modules=cythonize(extensions),
    include_dirs=[np.get_include()],
    cmdclass={'build_ext': build_ext_opt},
    python_requires='>=2.7',
    platforms='any'
)
