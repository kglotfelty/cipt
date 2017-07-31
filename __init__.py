from __future__ import absolute_import

"""
CIAO Image processing toolkit


This package provides a CIAOImage object that provides easy
access to the various CIAO image processing tools.

The underlying CIAO tools are run on temporary files (under the hood,
no user stuff).

"""

from .cipt import CIAOImage
from .enhanced_region import *

