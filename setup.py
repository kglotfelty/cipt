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
        url='https://github.com/kglotfelty/cipt/',
        packages=["ciao_contrib/cipt/"]
        )
"""
plot_shapes.py
setup.py
"""
