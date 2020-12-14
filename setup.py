#!/usr/bin/env python


import os
import sys

assert "ASCDS_INSTALL" in os.environ

from distutils.core import setup


setup( name='CIPT',
        version='4.13.0',
        description='CIAO Image Processing Toolkit',
        author='Kenny Glotfelty',
        author_email='glotfeltyk@si.edu',
        url='https://github.com/kglotfelty/cipy/',
        py_modules=["__init__",
                    "all",
                    "cipt",
                    "crateify",
                    "smooth_kernels",
                    "plot_shapes"]
        )
"""
plot_shapes.py
setup.py
"""
