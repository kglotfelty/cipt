#!/usr/bin/env python


import os
import sys

assert "ASCDS_INSTALL" in os.environ

ver = sys.version_info
os.environ["PYVER"] = "python{}.{}".format(ver[0],ver[1])

from distutils.core import setup


setup( name='CIPT',
        version='0.0.1',
        description='CIAO Image Processing Toolkit',
        author='Anonymous',
        author_email='kglotfelty@cfa.harvard.edu',
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
