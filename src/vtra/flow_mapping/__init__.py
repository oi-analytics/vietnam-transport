# -*- coding: utf-8 -*-
"""vtra
"""
import pkg_resources

__author__ = "Oxford Infrastructure Analytics and Contributors"
__copyright__ = "Oxford Infrastructure Analytics and Contributors"
__license__ = "mit"

try:
    __version__ = pkg_resources.get_distribution(__name__).version
except Exception:
    __version__ = 'unknown'
